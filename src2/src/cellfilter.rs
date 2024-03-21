use crate::prog_opts::GenPermitListOpts;
use crate::utils as afutils;
use ahash::{AHasher, RandomState};
use anyhow::{anyhow, Context};
use bstr::io::BufReadExt;
use itertools::Itertools;
use libradicl::exit_codes;
use libradicl::BarcodeLookupMap;
use needletail::bitkmer;
use num_format::{Locale, ToFormattedString};
use serde::Serialize;
use serde_json::json;
use slog::crit;
use slog::info;
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;
use indexmap::map::IndexMap;

#[derive(Debug, PartialEq, Eq)]
pub struct hit_info {
    chr: u32, // could be converted to u8, with a hashmap mapping u8 to chromosome name
    start: u32,
    end: u32,
    barcode: needletail::bitkmer::BitKmerSeq
    // rec_id: u64,
}

pub fn write_bed(hit_info_vec:&[hit_info],
                bed_path:&Path,
                chr_map:&IndexMap<String, u32>,
                hm:&HashMap<u64, u64, ahash::RandomState>,
                rev:&bool,
                bc_len:&usize            
            ) -> Result<(), Box<dyn std::error::Error>> {
            
        let mut writer = std::fs::File::create(bed_path)?;
        let keys: Vec<_> = chr_map.keys().cloned().collect();
        let chunk_size = 1000;
        let mut chunk = 0;
        let mut s = "".to_string() ;       
        for i in 0..hit_info_vec.len() {
            if chunk == chunk_size {
                chunk = 0;
            }
            let chr_val = keys[hit_info_vec[i].chr as usize].clone();
            let start = hit_info_vec[i].start;
            let end = hit_info_vec[i].end;
            let bc = hm.get(&hit_info_vec[i].barcode).unwrap();
            let bc_string = get_bc_string(bc, rev, bc_len);
            let s2 = [chr_val, start.to_string(), end.to_string(), bc_string].join("\t");
            s.push_str(&s2);
            s.push('\n');
    
            if chunk == chunk_size-1 {
                writer.write_all(s.as_bytes()).unwrap();
            }
            chunk+=1;
        }
        Ok(())       
}

impl Ord for hit_info {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        if self.barcode != other.barcode {
            self.barcode.cmp(&other.barcode)
        } else if self.start != other.start {
            self.start.cmp(&other.start)
        } else {
            self.end.cmp(&other.end)
        }
    }
}
impl PartialOrd for hit_info {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        if self.barcode != other.barcode {
            Some(self.barcode.cmp(&other.barcode))
        } else if self.start != other.start {
            Some(self.start.cmp(&other.start))
        } else {
            Some(self.end.cmp(&other.end))
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub enum CellFilterMethod {
    // use the distance method to
    // automatically find the knee
    // in the curve
    KneeFinding,
    // barcodes will be provided in the
    // form of an *unfiltered* external
    // permit list
    UnfilteredExternalList(PathBuf, usize),
}

pub fn get_bc_string(kmerseq: &needletail::bitkmer::BitKmerSeq,
    reverse_barcode: &bool,
    bc_len: &usize
) -> String {
    let kmseq = *kmerseq;
    let mut km: needletail::bitkmer::BitKmer = (kmseq, *bc_len as u8);
    if *reverse_barcode {
        km = needletail::bitkmer::reverse_complement(km);
    }
    let bytes = needletail::bitkmer::bitmer_to_bytes(km);
    let seq:String = String::from_utf8(bytes).expect("Invalid barcode");
    seq
}

pub fn update_barcode_hist_unfiltered(
    hist: &mut HashMap<u64, u64, ahash::RandomState>,
    unmatched_bc: &mut Vec<u64>,
    max_ambiguity_read: &mut usize,
    reverse_barcode: &bool,
    hit: &String,
    prev_hits: &mut u32,
    first_inst: &mut bool,
    hit_info_vec: &mut Vec<hit_info>,
    chr_map: &mut IndexMap<String, u32>,
    chr_count: &mut u32,
    uni_count: &mut u32,
    mult_count: &mut u32,
) -> usize {
    let mut new_read = 0usize;
    let split_line: Vec<String> = hit.split_whitespace().map(|s| s.to_string()).collect();
    let num_hits = split_line[4].parse::<u32>().unwrap();
    let bc = &split_line[3].as_bytes();
    let l = bc.len();
    // need a way to get rid of unwrap
    let mut km: needletail::bitkmer::BitKmer = needletail::bitkmer::BitNuclKmer::new(&bc[..], l as u8, false).next().unwrap().1;
    if *reverse_barcode {
        km = needletail::bitkmer::reverse_complement(km);
    }
    
    // Making sure we only count the multimapping hit once
    // If the number of hits in the prev record were 1, then current hit is the first one
    if *prev_hits == 1 {
        *first_inst = true;
        if num_hits > 1 {
            *mult_count += 1;
        }
        else {
            *uni_count += 1;
        }
    }
    // If the number of prev hits equal to current total, then
    // next hit will be the first one
    if *prev_hits == num_hits {
        *prev_hits = 1;
    } else {
        *prev_hits += 1;
    }

    let mut rec_id:i64 = -1;
    let mut prev_rec_id:i64 = -1;
    let nrecs = hit_info_vec.len();

    // if nrecs > 1 {
    //     prev_rec_id = hit_info_vec[nrecs-1].rec_id as i64;
    // }
    if *first_inst {
        if num_hits > 1 {
            *first_inst = false;
            *max_ambiguity_read = (num_hits as usize).max(*max_ambiguity_read);
        }
        // rec_id = prev_rec_id + 1; 
        // if let Some((_, km, _)) =
        //     needletail::bitkmer::BitNuclKmer::new(&bc[..], l as u8, false).next()
        // {
        //     let rc = needletail::bitkmer::reverse_complement(km);
        match hist.get_mut(&km.0) {
            // match hist.get_mut(&km.0) {
            Some(c) => *c += 1,
            None => {
                unmatched_bc.push(km.0);
            } // }
        }
        new_read = 1;
    } else {
        new_read = 0;
        rec_id = prev_rec_id;
    }
    let chr_s = split_line[0].clone();
    let mut chr = *chr_count;
    if chr_map.contains_key(&chr_s) {
        chr = *chr_map.get(&chr_s).unwrap();
    }
    else {
        chr_map.insert(chr_s.clone(), chr);
        *chr_count += 1;
    }

    if num_hits == 1{
        hit_info_vec.push(hit_info{
            chr: chr,
            start: split_line[1].parse::<u32>().unwrap(),
            end: split_line[2].parse::<u32>().unwrap(),
            barcode: km.0             
        });
    }

    new_read
}

fn populate_unfiltered_barcode_map<T: Read>(
    br: BufReader<T>,
    first_bclen: &mut usize,
) -> HashMap<u64, u64, ahash::RandomState> {
    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut hm = HashMap::with_hasher(s);

    // read through the external unfiltered barcode list
    // and generate a vector of encoded barcodes
    // let mut kv = Vec::<u64>::new();
    for l in br.byte_lines().flatten() {
        if *first_bclen == 0 {
            *first_bclen = l.len();
        } else {
            assert_eq!(
                *first_bclen,
                l.len(),
                "found barcodes of different lengths {} and {}",
                *first_bclen,
                l.len()
            );
        }
        if let Some((_, km, _)) =
            needletail::bitkmer::BitNuclKmer::new(&l[..], l.len() as u8, false).next()
        {
            // let rc = needletail::bitkmer::reverse_complement(km);
            hm.insert(km.0, 0);
            // hm.insert(rc.0, 0);
        }
    }
    hm
}

fn process_unfiltered(
    mut hm: &mut HashMap<u64, u64, ahash::RandomState>,
    mut unmatched_bc: Vec<u64>,
    bc_len: &usize,
    filter_meth: &CellFilterMethod,
    output_dir: &PathBuf,
    max_ambiguity_read: usize,
    cmdline: &str,
    log: &slog::Logger,
    gpl_opts: &GenPermitListOpts,
) -> anyhow::Result<u64> {
    let parent = std::path::Path::new(output_dir);
    std::fs::create_dir_all(parent)
        .with_context(|| format!("couldn't create directory path {}", parent.display()))?;

    // the smallest number of reads we'll allow per barcode
    let min_freq = match filter_meth {
        CellFilterMethod::UnfilteredExternalList(_, min_reads) => {
            info!(log, "minimum num reads for barcode pass = {}", *min_reads);
            *min_reads as u64
        }
        _ => {
            unimplemented!();
        }
    };

    // the set of barcodes we'll keep
    let mut kept_bc = Vec::<u64>::new();

    // iterate over the count map
    for (k, v) in hm.iter_mut() {
        // if this satisfies our requirement for the minimum count
        // then keep this barcode
        if *v >= min_freq {
            kept_bc.push(*k);
        } else {
            // otherwise, we have to add this barcode's
            // counts to our unmatched list
            for _ in 0..*v {
                unmatched_bc.push(*k);
            }
            // and then reset the counter for this barcode to 0
            *v = 0u64;
        }
    }

    // drop the absent barcodes from hm
    hm.retain(|_, &mut v| v > 0);

    // how many we will keep
    let num_passing = kept_bc.len();
    info!(
        log,
        "num_passing = {}",
        num_passing.to_formatted_string(&Locale::en)
    );

    // now, we create a second barcode map with just the barcodes
    // for cells we will keep / rescue.
    let bcmap2 = BarcodeLookupMap::new(kept_bc, *bc_len as u32);
    info!(
        log,
        "found {} cells with non-trivial number of reads by exact barcode match",
        bcmap2.barcodes.len().to_formatted_string(&Locale::en)
    );

    // finally, we'll go through the set of unmatched barcodes
    // and try to rescue those that have a *unique* neighbor in the
    // list of retained barcodes.

    //let mut found_exact = 0usize;
    let mut found_approx = 0usize;
    let mut ambig_approx = 0usize;
    let mut not_found = 0usize;

    let start_unmatched_time = Instant::now();

    unmatched_bc.sort_unstable();

    let mut distinct_unmatched_bc = 0usize;
    let mut distinct_recoverable_bc = 0usize;

    // mapping the uncorrected barcode to what it corrects to
    let mut corrected_list = Vec::<(u64, u64)>::with_capacity(1_000_000);

    for (count, ubc) in unmatched_bc.iter().dedup_with_count() {
        // try to find the unmatched barcode, but
        // look up to 1 edit away
        match bcmap2.find_neighbors(*ubc, false) {
            // if we have a match
            (Some(x), n) => {
                let cbc = bcmap2.barcodes[x];
                // if the uncorrected barcode had a
                // single, unique retained neighbor
                if cbc != *ubc && n == 1 {
                    // then increment the count of this
                    // barcode by 1 (because we'll correct to it)
                    if let Some(c) = hm.get_mut(&cbc) {
                        *c += count as u64;
                        corrected_list.push((*ubc, cbc));
                    }
                    // this counts as an approximate find
                    found_approx += count;
                    distinct_recoverable_bc += 1;
                }
                // if we had > 1 single-mismatch neighbor
                // then don't keep the barcode, but remember
                // the count of such events
                if n > 1 {
                    ambig_approx += count;
                }
            }
            // if we had no single-mismatch neighbor
            // then this barcode is not_found and gets
            // dropped.
            (None, _) => {
                not_found += count;
            }
        }
        distinct_unmatched_bc += 1;
    }
    let unmatched_duration = start_unmatched_time.elapsed();
    let num_corrected = distinct_recoverable_bc as u64;

    info!(
        log,
        "There were {} distinct unmatched barcodes, and {} that can be recovered",
        distinct_unmatched_bc,
        distinct_recoverable_bc
    );
    info!(
        log,
        "Matching unmatched barcodes to retained barcodes took {:?}", unmatched_duration
    );
    info!(log, "Of the unmatched barcodes\n============");
    info!(
        log,
        "\t{} had exactly 1 single-edit neighbor in the retained list",
        found_approx.to_formatted_string(&Locale::en)
    );
    info!(
        log,
        "\t{} had >1 single-edit neighbor in the retained list",
        ambig_approx.to_formatted_string(&Locale::en)
    );
    info!(
        log,
        "\t{} had no neighbor in the retained list",
        not_found.to_formatted_string(&Locale::en)
    );

    let parent = std::path::Path::new(output_dir);
    std::fs::create_dir_all(parent).with_context(|| {
        format!(
            "couldn't create path to output directory {}",
            parent.display()
        )
    })?;
    let o_path = parent.join("permit_freq.bin");

    match afutils::write_permit_list_freq(&o_path, *bc_len as u16, &hm) {
        Ok(_) => {}
        Err(error) => {
            panic!("Error: {}", error);
        }
    };

    /*
    // don't need this right now
    let s_path = parent.join("bcmap.bin");
    let s_file = std::fs::File::create(&s_path).expect("could not create serialization file.");
    let mut s_writer = BufWriter::new(&s_file);
    bincode::serialize_into(&mut s_writer, &bcmap2).expect("couldn't serialize barcode list.");
    */

    // now that we are done with using hm to count, we can repurpose it as
    // the correction map.
    for (k, v) in hm.iter_mut() {
        // each present barcode corrects to itself
        *v = *k;
    }
    for (uncorrected, corrected) in corrected_list.iter() {
        hm.insert(*uncorrected, *corrected);
    }

    let pm_path = parent.join("permit_map.bin");
    let pm_file = std::fs::File::create(pm_path).context("could not create serialization file.")?;
    let mut pm_writer = BufWriter::new(&pm_file);
    bincode::serialize_into(&mut pm_writer, &hm)
        .context("couldn't serialize permit list mapping.")?;

    let meta_info = json!({
    "max-ambig-record" : max_ambiguity_read,
    "cmd" : cmdline,
    "permit-list-type" : "unfiltered",
    "gpl_options" : &gpl_opts
    });

    let m_path = parent.join("generate_permit_list.json");
    let mut m_file = std::fs::File::create(m_path).context("could not create metadata file.")?;

    let meta_info_string =
        serde_json::to_string_pretty(&meta_info).context("could not format json.")?;
    m_file
        .write_all(meta_info_string.as_bytes())
        .context("cannot write to generate_permit_list.json file")?;

    info!(
        log,
        "total number of distinct corrected barcodes : {}",
        num_corrected.to_formatted_string(&Locale::en)
    );

    Ok(num_corrected)
}

pub fn generate_permit_list(gpl_opts: GenPermitListOpts) -> anyhow::Result<u64> {
    let bed_dir = gpl_opts.input_dir;
    let output_dir = gpl_opts.output_dir;
    let filter_meth = gpl_opts.fmeth.clone();
    let version = gpl_opts.version;
    let cmdline = gpl_opts.cmdline;
    let log = gpl_opts.log;

    let i_dir = std::path::Path::new(&bed_dir);

    // should we assume this condition was already checked
    // during parsing?
    if !i_dir.exists() {
        crit!(
            log,
            "the input bed path {} does not exist",
            bed_dir.display()
        );
        // std::process::exit(1);
        return Err(anyhow!("execution terminated unexpectedly"));
    }

    let mut first_bclen = 0usize;
    let mut unfiltered_bc_counts = None;
    let rc = true; // reverse complement barcode
    if let CellFilterMethod::UnfilteredExternalList(fname, _) = &filter_meth {
        let i_file = File::open(fname).context("could not open input file")?;
        let br = BufReader::new(i_file);
        unfiltered_bc_counts = Some(populate_unfiltered_barcode_map(br, &mut first_bclen));
        info!(
            log,
            "number of unfiltered bcs read = {}",
            unfiltered_bc_counts
                .as_ref()
                .unwrap()
                .len()
                .to_formatted_string(&Locale::en)
        );
    }

    let i_file = File::open(i_dir.join("map.bed")).context("could not open input bed file")?;
    let br = BufReader::new(i_file);
    let mut num_reads: usize = 0;

    // if dealing with the unfiltered type
    // the set of barcodes that are not an exact match for any known barcodes
    let mut unmatched_bc: Vec<u64>;
    // let mut num_orientation_compat_reads = 0usize;
    let mut max_ambiguity_read = 0usize;
    // Tracking if a unique or a multihit
    let mut first_inst = true;
    let mut hit_info_vec = Vec::<hit_info>::new();
    let mut ff = 0;
    let mut prev_hits: u32 = 1;
    let mut uni_count: u32 = 0;
    let mut multi_count: u32 = 0;
    let mut chr_count: u32 = 0;
    let mut chr_map:IndexMap<String, u32> = IndexMap::new();
    
    match filter_meth {
        CellFilterMethod::UnfilteredExternalList(_, _min_reads) => {
            unmatched_bc = Vec::with_capacity(10000000);
            // the unfiltered_bc_count map must be valid in this branch

            if let Some(mut hmu) = unfiltered_bc_counts {
                let start_unmatched_time = Instant::now();
                for l in br.lines() {
                    let l = l?;
                    ff += 1;
                    num_reads += update_barcode_hist_unfiltered(
                        &mut hmu,
                        &mut unmatched_bc,
                        &mut max_ambiguity_read,
                        &rc,
                        &l,
                        &mut prev_hits,
                        &mut first_inst,
                        &mut hit_info_vec,
                        &mut chr_map,
                        &mut chr_count,
                        &mut uni_count,
                        &mut multi_count
                    );
                    // num_reads += c.reads.len();
                }
                let unmatched_duration = start_unmatched_time.elapsed();
                println!("{:?}", unmatched_duration);
                println!("{:?}", hit_info_vec.len());
                println!("{:?}", hit_info_vec[10000]);
                let l = unmatched_bc.len();
                info!(
                    log,
                    "observed {} --- unmatched_bc {} --- max ambiguity read occurs in {} refs",
                    num_reads.to_formatted_string(&Locale::en),
                    l.to_formatted_string(&Locale::en),
                    // hdr.num_chunks.to_formatted_string(&Locale::en),
                    max_ambiguity_read.to_formatted_string(&Locale::en)
                );
                info!(
                    log,
                    "uni_count {} --- multi_count {} --- chr_count {}",
                    uni_count.to_formatted_string(&Locale::en),
                    multi_count.to_formatted_string(&Locale::en),
                    // hdr.num_chunks.to_formatted_string(&Locale::en),
                    chr_count.to_formatted_string(&Locale::en)
                );
                let corr = process_unfiltered(
                    &mut hmu,
                    unmatched_bc,
                    &first_bclen,
                    &filter_meth,
                    output_dir,
                    max_ambiguity_read,
                    cmdline,
                    log,
                    &gpl_opts,
                );
                info!(log, "sorting");
                let hh = &mut hit_info_vec[0..100];
                hh.sort_unstable();
                info!(log, "done sorting");
                let parent = std::path::Path::new(output_dir);
                let sorted_out_bed = parent.join("sorted_map.bed");
                let _=write_bed(&hh, &sorted_out_bed, &chr_map, &hmu, &rc, &first_bclen);
                corr
            } else {
                Ok(0)
            }
        }
        _ => Ok(0),
    }
}
