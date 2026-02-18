use std::fs::File;
use indicatif::{ProgressBar};
use std::collections::{BTreeMap, HashMap};

use num_format::{Locale, ToFormattedString};

use atty::is;

use std::io::{Write};
use std::time::{Duration, SystemTime};

use chrono::{DateTime, Utc};

use std::fmt;

/// MappingInfo captures all mapping data and is a way to easily copy this data over multiple analysis runs.
pub struct MappingInfo{
	/// reads that did not pass the filters
	pub quality:usize,
	pub length:usize,
	analyzed:usize,
	pub n_s:usize,
	pub poly_a:usize,
	/// reads that had no cell id
    pub no_sample:usize,
    /// reads that have no match in the geneIds object
    pub no_data:usize,
    /// reads with cell id and gene id
    pub ok_reads:usize,
    /// reads with cell_id - gene_id is not checked
    pub cellular_reads: usize,
    pub multimapper: usize,
    /// reads that are duplicates on the UMI level per cell and gene
    pub pcr_duplicates:usize,
    /// the amount of ok_reads after which to write a entry into the log file
   	pub split:usize,
   	/// the others are explained in the quantify_rhapsody.rs file.
    log_iter:usize,
    pub log_writer: Option<File>,
    pub min_quality:f32, 
    pub max_reads:usize, 
    pub local_dup:usize,
    pub total:usize,
    pub absolute_start: SystemTime,
    realtive_start: Option<SystemTime>,
    tmp_counter: Option<SystemTime>,
    pub single_processor_time: Duration,
    pub multi_processor_time: Duration,
    pub file_io_time: Duration,
    pub subprocess_time: Duration,
    pub reads_log: BTreeMap<String, usize >,
    pub error_counts: HashMap<String, usize>,  // To store error types and their counts
    // log should also print (if not likely to tty)
    std_out_is_tty: bool,
    pub hist:Vec<usize>,

}

impl fmt::Display for MappingInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // --- helpers ---
        fn pct(n: usize, d: usize) -> f32 {
            if d == 0 { 0.0 } else { (n as f32 / d as f32) * 100.0 }
        }
        fn dur_str(d: Duration) -> String {
            let (h, m, s, ms) = MappingInfo::split_duration(d);
            format!("{h} h {m} min {s} sec {ms} ms")
        }

        let analyzed = if self.analyzed == 0 { self.total } else { self.analyzed };
        let unknown = self.quality + self.length + self.n_s + self.poly_a;

        // Header
        writeln!(f, "MappingInfo")?;
        writeln!(f, "  started: {}", {
            // keep it cheap, avoid chrono if you want, but you already use it:
            let now: DateTime<Utc> = Utc::now();
            now.to_string()
        })?;
        writeln!(f)?;

        // Core counters
        writeln!(f, "Counts")?;
        writeln!(f, "  total reads        : {}", self.total.to_formatted_string(&Locale::en))?;
        writeln!(
            f,
            "  cellular reads     : {} ({:.2}% of total)",
            self.cellular_reads.to_formatted_string(&Locale::en),
            pct(self.cellular_reads, self.total)
        )?;
        writeln!(
            f,
            "  ok reads           : {} ({:.2}% of analyzed)",
            self.ok_reads.to_formatted_string(&Locale::en),
            pct(self.ok_reads, analyzed)
        )?;
        writeln!(
            f,
            "  no cell id         : {} ({:.2}% of total)",
            self.no_sample.to_formatted_string(&Locale::en),
            pct(self.no_sample, self.total)
        )?;
        writeln!(
            f,
            "  no gene id         : {} ({:.2}% of total)",
            self.no_data.to_formatted_string(&Locale::en),
            pct(self.no_data, self.total)
        )?;
        writeln!(
            f,
            "  multimapper        : {} ({:.2}% of analyzed)",
            self.multimapper.to_formatted_string(&Locale::en),
            pct(self.multimapper, analyzed)
        )?;
        writeln!(
            f,
            "  pcr duplicates     : {}",
            self.pcr_duplicates.to_formatted_string(&Locale::en)
        )?;
        writeln!(
            f,
            "  filtered (unknown) : {} ({:.2}% of total)",
            unknown.to_formatted_string(&Locale::en),
            pct(unknown, self.total)
        )?;
        writeln!(
            f,
            "    -> bad quality   : {} ({:.2}% of total)",
            self.quality.to_formatted_string(&Locale::en),
            pct(self.quality, self.total)
        )?;
        writeln!(
            f,
            "    -> too short     : {} ({:.2}% of total)",
            self.length.to_formatted_string(&Locale::en),
            pct(self.length, self.total)
        )?;
        writeln!(
            f,
            "    -> N's           : {} ({:.2}% of total)",
            self.n_s.to_formatted_string(&Locale::en),
            pct(self.n_s, self.total)
        )?;
        writeln!(
            f,
            "    -> polyA         : {} ({:.2}% of total)",
            self.poly_a.to_formatted_string(&Locale::en),
            pct(self.poly_a, self.total)
        )?;
        writeln!(f)?;

        // Read-type breakdown (if any)
        if !self.reads_log.is_empty() {
            writeln!(f, "Read types")?;
            for (name, value) in &self.reads_log {
                let denom = self.cellular_reads;
                writeln!(
                    f,
                    "  {:<20} {} reads ({:.2}% of cellular)",
                    format!("{name}:"),
                    value.to_formatted_string(&Locale::en),
                    pct(*value, denom)
                )?;
            }
            writeln!(f)?;
        }

        // Error report (if any)
        if !self.error_counts.is_empty() {
            writeln!(f, "Reported issues")?;
            writeln!(f, "  {:<32} {}", "Error Type", "Count")?;
            writeln!(f, "  {}", "-".repeat(32 + 1 + 12))?;
            // stable order: sort by key
            let mut keys: Vec<_> = self.error_counts.keys().collect();
            keys.sort();
            for k in keys {
                let c = self.error_counts.get(k).copied().unwrap_or(0);
                writeln!(
                    f,
                    "  {:<32} {}",
                    k,
                    c.to_formatted_string(&Locale::en)
                )?;
            }
            writeln!(f)?;
        }

        // Histogram (only show if anything non-zero)
        if self.hist.iter().any(|&x| x != 0) {
            writeln!(f, "Histogram")?;
            for (i, &v) in self.hist.iter().enumerate() {
                if v != 0 {
                    writeln!(f, "  bin {:>2}: {}", i, v.to_formatted_string(&Locale::en))?;
                }
            }
            writeln!(f)?;
        }

        // Timings
        writeln!(f, "Timings")?;
        writeln!(f, "  overall     : {}", dur_str(self.absolute_start.elapsed().unwrap_or(Duration::new(0, 0))))?;
        writeln!(f, "  file I/O    : {}", dur_str(self.file_io_time))?;
        writeln!(f, "  single-cpu  : {}", dur_str(self.single_processor_time))?;
        writeln!(f, "  multi-cpu   : {}", dur_str(self.multi_processor_time))?;
        if self.subprocess_time != Duration::new(0, 0) {
            writeln!(f, "  subprocess  : {}", dur_str(self.subprocess_time))?;
        }

        Ok(())
    }
}

impl MappingInfo{
	pub fn new(log_writer:Option<File>, min_quality:f32, max_reads:usize, ) -> Self{
		let absolute_start = SystemTime::now();
		let single_processor_time = Duration::new(0,0);
		let multi_processor_time = Duration::new(0,0);
		let file_io_time = Duration::new(0,0);
		let reads_log = BTreeMap::new();
		let subprocess_time = Duration::new(0,0);
		let mut this = Self{
			quality: 0,
		    length: 0,
		    analyzed: 1,
		    n_s: 0,
		    poly_a: 0,
			no_sample: 0,
			no_data: 0,
			ok_reads: 0,
			cellular_reads: 0,
			multimapper: 0,
			pcr_duplicates: 0,
			split: 1_000_000,
			log_iter: 0,
			log_writer,
			min_quality,
			max_reads,
			local_dup: 0,
			total: 0,
			absolute_start,
			realtive_start: None,
			tmp_counter: None,
			single_processor_time,
			multi_processor_time,
			file_io_time,
			subprocess_time,
			reads_log,
			error_counts: HashMap::new(),  // Initialize the HashMap
			std_out_is_tty: is(atty::Stream::Stdout) ,
			hist: vec![0; 20],
		};
		this.start_counter();
		this
	}

	pub fn iterate_hist( &mut self, id: usize) {
		if id < self.hist.len(){
			self.hist[id] +=1;
		}
	}

	pub fn start_counter ( &mut self ){
		self.realtive_start = Some( SystemTime::now() );
	}

	pub fn start_ticker ( &mut self )  {
		self.tmp_counter = Some( SystemTime::now() );
	}
	pub fn stop_ticker ( &mut self ) -> ( u128, u128, u128, u128 ) {
		let ret = MappingInfo::split_duration( self.tmp_counter.unwrap_or( SystemTime::now() ).elapsed().unwrap() );
		self.tmp_counter = Some( SystemTime::now() );
		ret
	}

	pub fn stop_single_processor_time ( &mut self ) {
		self.single_processor_time += self.realtive_start.unwrap().elapsed().unwrap();
		self.start_counter();
	}

	pub fn stop_multi_processor_time ( &mut self ) {
		self.multi_processor_time += self.realtive_start.unwrap().elapsed().unwrap();
		self.start_counter();
	}

	pub fn subprocess_time ( &mut self ) {
		self.subprocess_time += self.realtive_start.unwrap().elapsed().unwrap();
		self.start_counter();
	}

	pub fn stop_file_io_time ( &mut self ) {
		self.file_io_time += self.realtive_start.unwrap().elapsed().unwrap();
		self.start_counter();
	}

	pub fn elapsed_time_split ( &self ) -> ( u128, u128, u128, u128 ){
		let elapsed = self.absolute_start.elapsed().unwrap();
		MappingInfo::split_duration( elapsed )
	}

	pub fn now (&self) -> String{

		let now: DateTime<Utc> = Utc::now();
		format!("{}", now)
	    
	}

	// Unified reporting method that logs errors into the HashMap
    pub fn report<S: AsRef<str>>(&mut self, issue: S) {
        let issue = issue.as_ref();
        // Increment the count for the issue type
        *self.error_counts.entry(issue.to_string()).or_insert(0) += 1;
        //println!("Issue reported: {}", issue);
    }

    // Optionally, add a method to retrieve counts for a specific issue
    pub fn get_issue_count(&self, issue: &str) -> usize {
        *self.error_counts.get(issue).unwrap_or(&0)  // Return count or 0 if not present
    }

    // Method to export error_counts to a CSV file
    pub fn report_to_csv(&self, file_path: &str) {
        let mut file = File::create(file_path).unwrap();  // Create a file for writing
        writeln!(file, "Report Type\tCount").unwrap();  // Write CSV header

        // Iterate over the error_counts and write each as a row in the CSV
        for (error_type, count) in &self.error_counts {
            match writeln!(file, "{}\t{}", error_type, count){
            	Ok(_) => {},
            	Err(err) => eprintln!("An error occured while exporting the mapping report:\n{err}"),
            };  // Write each error type and count
        }

        //Ok(())  // Return Ok if successful
    }

    // Method to export error_counts to a CSV-formatted String
	pub fn report_to_string(&self) -> String {
	    // Start with the header
	    let mut output = String::from("Error Type\tCount\n");

	    // Iterate over the error_counts and append each as a row in the CSV format
	    for (error_type, count) in &self.error_counts {
	        // Append each error type and count, followed by a newline
	        let formatted = count.to_formatted_string(&Locale::en);
	        output.push_str(&format!("{}\t{}\n", error_type, formatted));
	    }

	    output // Return the content as a String
	}

	pub fn split_duration( elapsed:Duration ) -> ( u128, u128, u128, u128 ){

        let mut milli = elapsed.as_millis();

        let mil = milli % 1000;
        milli= (milli - mil) /1000;

        let sec = milli % 60;
        milli= (milli -sec) /60;

        let min = milli % 60;
        milli= (milli -min) /60;

        (milli, min, sec, mil )

    }

    pub fn iter_read_type(&mut self, name:&str ){
    	*self.reads_log.entry(name.to_string()).or_insert(0) += 1; 
    }

    pub fn read_types_to_string(&self, names:Vec<&str> ) -> String {
        let mut formatted_entries = String::new();

        for name in &names {
            let value = self.reads_log.get(*name).unwrap_or(&0);
        	let formatted_name = format!("{:<18}:", name); 
            formatted_entries.push_str(&format!("{} {} reads ({:.2}% of cellular)\n", formatted_name, value, *value as f32 / self.cellular_reads as f32 *100_f32 ));
        }
        formatted_entries
    }


	pub fn merge(&mut self, other:&MappingInfo ){
		self.no_sample += other.no_sample;
		self.no_data += other.no_data;
		//unknown is defined without multiprocessor support
		self.ok_reads += other.ok_reads;
		self.pcr_duplicates += other.pcr_duplicates;
		self.cellular_reads += other.cellular_reads;
		for (name, value) in &other.reads_log {
			*self.reads_log.entry(name.to_string()).or_insert(0) += value;
		}
		self.analyzed = self.total;
		for (error_type, count) in &other.error_counts {
            // For each error type in `other`, increment the value in `self`
            *self.error_counts.entry(error_type.clone()).or_insert(0) += count;
        }
        for (a, b) in self.hist.iter_mut().zip(&other.hist) {
		    *a += *b;
		}
	}



	pub fn write_to_log ( &mut self, text:String ){

		match &mut self.log_writer{
			Some(file) => {
				match writeln!( file , "{text}" ){
		            Ok(_) => (),
		            Err(err) => {
		                eprintln!("write error: {err}" );
		            }
		        };
			},
			None => {},
		}
		
	}


	/// add info from .report() into the log file as tab sep tables.
	pub fn log_report( &mut self ) {
		let log_str = self.report_to_string();
		self.write_to_log( log_str );
	}

	pub fn log( &mut self, pb:&ProgressBar ){
		if self.total % self.split == 0{
			self.log_iter+=1;
            let log_str = self.log_str();
            pb.set_message( log_str.clone() );
            pb.inc(1);
            self.write_to_log( log_str );
            self.local_dup = 0;
		}
	}

	pub fn log_str( &mut self ) -> String{
		format!("{:.2} mio reads ({:.2}% with cell_id, {:.2}% with gene_id {:.2}% multimapper)",
            self.total as f32 / self.split as f32,
            self.cellular_reads as f32 / (self.analyzed) as f32 * 100.0 , 
            self.ok_reads as f32 / (self.analyzed) as f32 * 100.0,
            self.multimapper as f32 / (self.analyzed) as f32 * 100.0,
         )
	}
	pub fn program_states_string( &self ) -> String{
		let mut result = String::from("");
		let  (mut hours,mut min,mut sec ,mut mulli ) = Self::split_duration( self.absolute_start.elapsed().unwrap() );
	   	result += format!("   overall run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	( hours, min, sec , mulli ) = Self::split_duration( self.file_io_time);
	   	result += format!("   file-io run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	( hours, min, sec , mulli ) = Self::split_duration( self.single_processor_time);
	   	result += format!("single-cpu run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	( hours, min, sec , mulli ) = Self::split_duration( self.multi_processor_time);
	   	result += format!(" multi-cpu run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	if self.subprocess_time != Duration::new(0,0) {
	   		( hours, min, sec , mulli ) = Self::split_duration( self.subprocess_time);
	    	result += format!("subprocess run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	}
	   	result
	}

	pub fn summary( &mut self, reads_genes:usize, reads_ab :usize, reads_samples:usize ) -> String{

		let pcr_duplicates = self.cellular_reads - reads_genes - reads_ab - reads_samples;

		let unknown = self.quality + self.length + self.n_s + self.poly_a;
	    let mut result = "\nSummary:\n".to_owned()
	    	+format!(     "cellular   reads  : {} reads ({:.2}% of total)\n", self.cellular_reads, (self.cellular_reads as f32 / self.total as f32) * 100.0 ).as_str()
	    	+format!(     "no cell ID reads  : {} reads ({:.2}% of total)\n", self.no_sample, (self.no_sample as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     "no gene ID reads  : {} reads ({:.2}% of total)\n", self.no_data.saturating_sub(self.no_sample), ( self.no_data.saturating_sub( self.no_sample) as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     "filtered   reads  : {} reads ({:.2}% of total)\n", unknown, (unknown as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     " ->  multimapper  : {} reads ({:.2}% of total)\n", self.multimapper, ( self.multimapper as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     " -> bad qualiity  : {} reads ({:.2}% of total)\n", self.quality, ( self.quality as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     " ->    too short  : {} reads ({:.2}% of total)\n", self.length, ( self.length as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     " ->          N's  : {} reads ({:.2}% of total)\n", self.n_s, ( self.n_s as f32 / self.total as f32) * 100.0).as_str()
	    	+"\n"
	    	+format!(     "total      reads  : {} reads\n", self.total ).as_str()
	    	+"\ncollected read counts:\n"
	    	+self.read_types_to_string(vec!["expression reads", "antibody reads", "sample reads"]).as_str()
	    	+"\nreported UMI counts:\n"
	    	+format!(     "expression reads  : {} UMIs ({:.2}% of cellular)\n", reads_genes, (reads_genes as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "antibody reads    : {} UMIs ({:.2}% of cellular)\n", reads_ab, (reads_ab as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "sample reads      : {} UMIs ({:.2}% of cellular)\n", reads_samples, (reads_samples as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	//+format!(     "unique reads      : {} reads ({:.2}% of cellular)\n", reads_genes + reads_ab + reads_samples, ( (reads_genes + reads_ab + reads_samples) as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "\nPCR duplicates or bad cells: {} reads ({:.2}% of cellular)\n\n", pcr_duplicates, ( pcr_duplicates as f32 / self.cellular_reads as f32 ) * 100.0 ).as_str()
	   		+"timings:\n";
	   	result += &self.program_states_string();
	   	self.write_to_log( result.clone() );
        result
	}
	
}
