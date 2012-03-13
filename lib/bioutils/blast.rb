##
# A set of methods focused on Blast (formatting, local/remote blast etc) that can be included via include Blast
# 
# @author Anthony Underwood
module Blast
  require 'utils/hash_reverse_merge'
  require 'utils/method_argument_parser'
  require 'tempfile'

  require 'bioutils/rich_sequence_utils'
  require 'bioutils/ncbi_blast'


  require 'rubygems'
  require 'bio'
  require 'json'
  # A method to produce a local blast database based on a directory containing genbank/embl files
  # @param [Hash] options The options for the method
  # @option options [String] :output_dir Path to the output directory (default is .)
  # @option options [String] :formatdb_options The command line arguments that wil be used with formatdb
  # @option options [String] :input_dir The path to the directory containing the rich sequence files including the glob e.g *.gbk
  # @option options [String] :database_name The name that will be used for the prefix for the blast database
  # @option options [String] :database_labels The type of label that will be used to create
  # the fasta files used to create a database, either :counter or :full_annotation,
  # :simple_gene_and_product, :gene_and_product
  # @option options [String] :final_db_location The path to the final db location
  def blast_database_from_directory(options)
    default_options = {
      :output_dir => ".",
      :formatdb_options => "-o T -p F",
      :formatdb_dir => "/usr/local/blast/bin/",
      :database_labels => :counter,
      :final_db_location => "."
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :input_dir, :required => true, :type => :string
    end

    tmpfile_object = Tempfile.new('temp')
    
    rich_sequence_files  = Dir.glob(options[:input_dir])
    counter = 0
    rich_sequence_files.each do |rich_sequence_file|
      counter += 1
      puts counter.to_s+" "+rich_sequence_file
      rich_sequence_object = rich_sequence_object_from_file(rich_sequence_file)
      write_blast_db_compatible_fasta_file(rich_sequence_object, tmpfile_object, options[:database_labels],options[:formatdb_options] =~ /-p T/)
    end
    tmpfile_object.close
    Dir.chdir(options[:output_dir])
    formatdb(:path_to_fasta_input_sequence => tmpfile_object.path, :formatdb_options => options[:formatdb_options], :formatdb_dir => options[:formatdb_dir], :database_name => options[:database_name],:final_db_location  => options[:final_db_location])
    return options[:database_name]
  end

  # A method to produce a local blast database
  # @param [Hash] options The options for the method
  # @option options [String] :output_dir Path to the output directory (default is .)
  # @option options [String] :formatdb_options The command line arguments that wil be used with formatdb
  # @option options [Array] :input_sequences An array of either paths to rich sequence files or Bioruby Bio::Sequence objects (can be mixed)
  # @option options [String] :database_name The name that will be used for the prefix for the blast database
  # @option options [String] :database_labels The type of label that will be used to create
  # the fasta files used to create a database, either :counter or :full_annotation,
  # :simple_gene_and_product, :gene_and_product
  # @option options [String] :final_db_location The path to the final db location
  def blast_database_from_rich_sequences(options)
    default_options = {
      :output_dir => ".",
      :formatdb_options => "-o T -p F",
      :formatdb_dir => "/usr/local/blast/bin",
      :database_labels => :counter,
      :final_db_location => ".",
      :protein => false
    }
    options[:formatdb_options] = "-o T -p T" if options[:protein]
    options.reverse_merge!(default_options)


    options = MethodArgumentParser::Parser.check_options options  do
      option :input_sequences, :required => true, :type => :array
    end
    # convert input_sequences to array if necessary
    options[:input_sequences] = Array(options[:input_sequences])
    # default database_name
    raise ArgumentError, "You must supply a database name (options[:database_name]) if you are using more than 1 sequence to format the database" if options[:input_sequences].size > 1 && options[:database_name].nil?
    options[:database_name] = File.basename(options[:input_sequences].first).sub(/#{File.extname(options[:input_sequences].first)}$/,'') if options[:database_name].nil?

    rich_sequence_objects = Array.new
    options[:input_sequences].each do |input_sequence|
      if input_sequence.class.to_s == "String"
        rich_sequence_object = rich_sequence_object_from_file(input_sequence)
      elsif input_sequence.class.to_s =~ "Bio::(Genbank|EMBL)"
        rich_sequence_object = input_sequence
      end
      if rich_sequence_object.is_a?(Array) # multiple rich sequence objects have been returned e.g a multi genbank file
        rich_sequence_objects += rich_sequence_object
      else
        rich_sequence_objects << rich_sequence_object
      end
    end
    blast_database_from_biosequences(:rich_sequence_objects => rich_sequence_objects,:formatdb_options => options[:formatdb_options], :formatdb_dir => options[:formatdb_dir],  :database_name => options[:database_name], :database_labels => options[:database_labels], :final_db_location  => options[:final_db_location], :protein => options[:protein])
  end

  # Create a blast database from a rich sequence (e.g genbank or embl) that has been used to make a Bio::Sequence object. N.B Need to make sure formatdb is in the PATH
  # @param [Hash] options The options for the method
  # @option options [Array] :rich_sequence_objects An array of either Bio::(Genbank|EMBL) objects that will be used to generate the database
  # @option options [String] :output_dir The path to where the database will be created
  # @option options [String] :formatdb_options The command line arguments that wil be used with formatdb
  # @option options [String] :database_name The prefix that wil be used to 'label' the database files
  # @option options [String] :database_labels The type of label that will be used to create the fasta files used to create a database, either :counter or :full_annotation
  # @option options [String] :final_db_location The path to the final db location
  # @return [String] The name of the blast database created
  def blast_database_from_biosequences(options)
    default_options = {
      :output_dir => ".",
      :formatdb_options => "-o T -p F",
      :formatdb_dir => "/usr/local/blast/bin",
      :database_labels => :counter,
      :final_db_location => "."
    }
    options.reverse_merge!(default_options)

    # create temp file for fasta file that will be used to create query blast database
    tmpfile_object = Tempfile.new('temp')
    options[:rich_sequence_objects].each do |rich_sequence_object|
      write_blast_db_compatible_fasta_file(rich_sequence_object, tmpfile_object, options[:database_labels], options[:protein])
    end
    tmpfile_object.close
    Dir.chdir(options[:output_dir])
    formatdb(:path_to_fasta_input_sequence => tmpfile_object.path, :formatdb_options => options[:formatdb_options], :formatdb_dir => options[:formatdb_dir], :database_name => options[:database_name],:final_db_location  => options[:final_db_location])
    return options[:database_name]
  end
  # a method to write a blast database compatible fasta file from a rich sequence object
  # @param [Bio::Sequence] rich_sequence_object A BioRuby rich sequence biosequence object
  # @param [File] file_object A ruby file handle
  def write_blast_db_compatible_fasta_file(rich_sequence_object, file_object, database_labels, protein = false)
    sequence = rich_sequence_object.seq
    feature_counter = 0
    rich_sequence_object.each_cds do |feature|
      feature_counter += 1
      location_string =  feature.locations.to_s
      predicted_feature_sequence = sequence.splice(location_string)
      predicted_feature_sequence = predicted_feature_sequence.translate(1,11) if protein
      database_id = create_local_blast_database_id(rich_sequence_object)
      if database_labels == :counter
        file_object.puts ">lcl|#{database_id}-#{feature_counter}" # use the lcl|id format as specified http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/formatdb_fastacmd.html/http://bioinformatics.abc.hu/tothg/biocomp/other/blast/formatdb.html
      elsif database_labels == :simple_gene_and_product
        qualifiers_array = feature.qualifiers.map do |qualifier|
          [qualifier.qualifier, qualifier.value] if ["product", "gene, locus_tag"].include?(qualifier.qualifier)
        end
        qualifiers_array.compact!
        qualifiers_array.unshift(["organism", "#{rich_sequence_object.definition} (#{rich_sequence_object.entry_id})"])
        annotation_simple = qualifiers_array.map {|pair| "#{pair[0]}: #{pair[1]}"}.join(", ")
        file_object.puts ">lcl|#{database_id}-#{feature_counter}  #{annotation_simple}" # add annotation as a JSON string
      elsif database_labels == :gene_and_product
        qualifiers_array = feature.qualifiers.map do |qualifier|
          [qualifier.qualifier, qualifier.value] if ["product", "gene, locus_tag"].include?(qualifier.qualifier)
        end
        qualifiers_array.compact!
        qualifiers_array.unshift(["location", location_string])
        qualifiers_array.unshift(["organism", "#{rich_sequence_object.definition} (#{rich_sequence_object.entry_id})"])
        annotation_json = qualifiers_array.to_json
        file_object.puts ">lcl|#{database_id}-#{feature_counter}  #{annotation_json}" # add annotation as a JSON string
      elsif database_labels == :full_annotation
        qualifiers_array = feature.qualifiers.map{|qualifier| [qualifier.qualifier, qualifier.value]}
        qualifiers_array.unshift(["location", location_string])
        qualifiers_array.unshift(["organism", "#{rich_sequence_object.definition} (#{rich_sequence_object.entry_id})"])
        annotation_json = qualifiers_array.to_json
        file_object.puts ">lcl|#{database_id}-#{feature_counter}  #{annotation_json}" # add annotation as a JSON string
      end
      file_object.puts predicted_feature_sequence
    end
  end
  # a method to format a blast database
  # @param [Hash] options The options for the method
  # @option options [String] :path_to_fasta_input_sequence The path to the multifasta file that will be used to format the database
  # @option options [String] :formatdb_dir The directory in which to find formatdb MUST include trailing slash (to be used if formatdb not in the path)
  # @option options [String] :formatdb_options The command line arguments that wil be used with formatdb
  # @option options [String] :database_name The prefix that wil be used to 'label' the database files
  # @option options [String] final_db_location The path to the final db location
  def formatdb(options)
    default_options = {
      :formatdb_options => "-o T -p F",
      :formatdb_dir => "/usr/local/blast/bin",
      :final_db_location => "."
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :path_to_fasta_input_sequence, :required => true, :type => :string
      option :database_name,  :required => true, :type => :string
    end
    `#{options[:formatdb_dir]}/formatdb -i "#{options[:path_to_fasta_input_sequence]}" #{options[:formatdb_options]} -n #{options[:database_name]}`
    unless options[:final_db_location]  == "."
      require 'fileutils'
      Dir.mkdir(options[:final_db_location]) unless File.exists?(options[:final_db_location])
      if options[:formatdb_options] =~ /-p T/
        FileUtils.move Dir.glob("#{options[:database_name]}*.p*"), "#{options[:final_db_location]}/"
      else
        FileUtils.move Dir.glob("#{options[:database_name]}*.n*"), "#{options[:final_db_location]}/"
      end
    end
  end

  # a method to format a database based on supplied input sequence(s)
  # @param [Hash] options The options for the method
  # @option options [Array] input_sequences An array of strings which are the path(s) to the input sequences
  # @option options [String] :formatdb_dir The directory in which to find formatdb MUST include trailing slash (to be used if formatdb not in the path)
  # @option options [String] :formatdb_options The command line arguments that wil be used with formatdb
  # @option options [String] :database_name The prefix that wil be used to 'label' the database files
  # @option options [String] :final_db_location The path to the final db location
  def blast_database_from_sequences(options)
    default_options = {
      :formatdb_options => "-o T -p F",
      :formatdb_dir => "/usr/local/blast/bin/",
      :final_db_location => "."
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :input_sequences, :required => true, :type => :array
      option :database_name,  :required => true, :type => :string
    end
    fasta_sequence_paths = Array.new
    temp_files = Array.new
    options[:input_sequences].each do |input_sequence|
      input_sequence_format = guess_sequence_format(input_sequence)
      if input_sequence_format != :fasta
        temp_files << Tempfile.new('temp')
        rich_sequence_to_fasta(input_sequence, temp_files.last)
        temp_files.last.close
        fasta_sequence_paths << temp_files.last.path
      else
        fasta_sequence_paths << input_sequence
      end
    end

    formatdb(:path_to_fasta_input_sequence => fasta_sequence_paths.join(" "), :formatdb_options => options[:formatdb_options], :formatdb_dir => options[:formatdb_dir], :database_name => options[:database_name], :final_db_location  => options[:final_db_location])
  end

  # A method to perform a local blast search
  # @param [Hash] options The options for the method
  # @option options [String,Bio::Sequence] :sequence The query sequence either as a string or as a Bio::Sequence object
  # @option options [String] :blast_program The blast program (blastn is the default)
  # @option options [String] :blast_database The path and name of the local blast database
  # @option options [String] :path_to_blast_executable The path to the blast executables
  # @option options [String] :blast_options Blast options e.g default is '-e 1e-10' (see http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastall/blastall_node23.html)
  def blast_sequence_locally(options)
    default_options = {
      :blast_program => "blastn",
      :path_to_blast_executable => "/usr/local/blast/bin/blastall",
      :blast_options => "-e 1e-10 -F F"
    }
    options.reverse_merge!(default_options)
    options = MethodArgumentParser::Parser.check_options options  do
      option :blast_database, :required => true, :type => String
    end
    
    raise ArgumentError, "blastall can not be found at #{options[:path_to_blast_executable]}" unless File.exists?(options[:path_to_blast_executable])
    local_blast_factory = Bio::Blast.local(options[:blast_program], options[:blast_database],options[:blast_options],options[:path_to_blast_executable])
    attempts = 0
    begin
      report = local_blast_factory.query(options[:sequence])
    rescue
      puts "local blast caused an error (blast program: #{options[:blast_program]}, blast database: #{options[:blast_database]}, blast options: #{options[:blast_options]}, path to blast executable: #{options[:path_to_blast_executable]}"
      attempts += 1
      if attempts < 6 #5 retries
        retry
      else
        return nil
      end
    end
    return report
  end
    # A method to perform a local blast search without using bioruby blast factory
  # @param [Hash] options The options for the method
  # @option options [String,Bio::Sequence] :sequence The query sequence either as a string or as a Bio::Sequence object
  # @option options [String] :blast_program The blast program (blastn is the default)
  # @option options [String] :blast_database The path and name of the local blast database
  # @option options [String] :path_to_blast_executable The path to the blast executables
  # @option options [String] :blast_options Blast options e.g default is '-e 1e-10' (see http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastall/blastall_node23.html)
  def blast_sequence_locally_without_bioruby(options)
    default_options = {
      :blast_program => "blastn",
      :path_to_blast_executable => "/usr/local/blast/bin/blastall",
      :blast_options => "-e 1e-10 -F F"
    }
    options.reverse_merge!(default_options)
    options = MethodArgumentParser::Parser.check_options options  do
      option :blast_database, :required => true, :type => String
    end
    
    
    raise ArgumentError, "blastall can not be found at #{options[:path_to_blast_executable]}" unless File.exists?(options[:path_to_blast_executable])
    
    tmpfile_object = Tempfile.new('temp')
    tmpfile_object.puts options[:sequence].seq
    tmpfile_object.close
    attempts = 0
    report = nil

    begin
      report_text = ""
      IO.popen("#{options[:path_to_blast_executable]} #{options[:blast_options]} -p #{options[:blast_program]} -d #{options[:blast_database]} -i #{tmpfile_object.path} -m 7 2>&1") do |cmd|
        cmd.each do |line|
          if line =~ /error for object/
            Process.kill 'TERM', cmd.pid
            report_text = nil
            break
          else
            report_text += line
          end
        end
        cmd.close
        report = Bio::Blast::Report.new(report_text) unless report_text.nil?
      end
    rescue
      puts "local blast caused an error #{$! }(blast command: #{options[:path_to_blast_executable]} #{options[:blast_options]} -p #{options[:blast_program]} -d #{options[:blast_database]} -i #{tmpfile_object.path} -m 7 2>&1, current directory: #{Dir.pwd})"
      attempts += 1
      if attempts < 6 #5 retries
        retry
      else
        return nil
      end
    end
    tmpfile_object.delete
    return report
  end
  # A method to perform a remote blast search
  # @param [Hash] options The options for the method
  # @option options [String,Bio::Sequence] :sequence The query sequence either as a string or as a Bio::Sequence object
  # @option options [String] :blast_program The blast program (blastn is the default)
  # @option options [String] :blast_database The remote database to be searched ('nr' by default)
  # @option options [Hash] :blast_options A hash of options to fine tune the blast search e.g :expect => '1e-3', :filter => 'T',
  # @option options [Hash] :retrieve_options A hash of options to fine tune the blast result retrieval. defaults are 'DESCRIPTIONS' => 100,'FORMAT_TYPE' => 'Text', 'ALIGNMENTS' => 50. Can change FORMAT TYPE to XML
  def blast_sequence_remotely(options)
    default_options = {
      :blast_program => "blastn",
      :blast_database => "nr",
      :blast_options => {},
      :retrieve_options => {}
    }
    options.reverse_merge!(default_options)

    remote_blast_factory = Bio::Blast::NCBI.new(:blast_program => options[:blast_program], :blast_database => options[:blast_database], :submit_params => options[:blast_options], :retrieve_params => options[:retrieve_options])
    attempts = 0
    begin
      report = remote_blast_factory.query(options[:sequence])
    rescue
      attempts += 1
      if attempts < 6 #5 retries
        puts "remote blast error"
        retry
      else
        return nil
      end
    end
    return report
  end
  # A method to determine if a hit is significant or not
  # @param [Hash] options The options for the method
  # @option options [Bio::Blast::Report::Hit] :hit The hit to be tested
  # @option options [Float] :evalue_cutoff The cutoff to determine whether a hit is significant or not
  # @option options [Float] :minimum_hit_length The minumum percentage overlap between the query and the subject hit to determine if it is significant
  # @option options [Float] :percent_identity_cutoff The minumum percentage identity between the query and the subject hit to determine if it is significant
  # @option options [String,Bio::Sequence] :query_sequence The query sequence either as a string or as a Bio::Sequence object
  # @return [Boolean] Returns either true or false
  def significant_hit?(options)
    default_options = {
      :percent_identity_cutoff => 90,
      :minimum_hit_length => 85,
      :display_metrics => false
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :hit, :required => true
      option :query_sequence,  :required => true
    end
    
    
    # calculate total hsp length
    total_length_of_hsps, total_number_of_matches = determine_hsp_length_and_mathes(options[:hit])
    
    hit_length = options[:hit].len
    hit_length = options[:hit].hsps.first.align_len if hit_length.nil? # if a remote blast hit them len will not be defined

    
    shortest_sequence = [options[:query_sequence].length,hit_length].min
    total_length_of_hsps_as_percentage_of_shortest = (total_length_of_hsps/shortest_sequence.to_f)*100

    if options[:display_metrics]
      puts "percent identity: #{total_number_of_matches/total_length_of_hsps.to_f*100 }. percentage length matched: #{total_length_of_hsps_as_percentage_of_shortest}. evalue: #{options[:hit].evalue}"
    end

    if options[:evalue_cutoff] # by default based on percent id but if evalue specified will use this as a cutofff
      if options[:hit].evalue < options[:evalue_cutoff] && total_length_of_hsps_as_percentage_of_shortest >= options[:minimum_hit_length]# sig hit
        return true
      else
        return false
      end
    else
      if total_number_of_matches/total_length_of_hsps.to_f*100 >= options[:percent_identity_cutoff] && total_length_of_hsps_as_percentage_of_shortest >= options[:minimum_hit_length]# sig hit
        return true
      else
        return false
      end
    end
  end

  def determine_hsp_length_and_mathes(hit)
    total_length_of_hsps = 0
    total_number_of_matches = 0

    hit.hsps.each do |hsp|
      total_length_of_hsps += (hsp.query_to - hsp.query_from).abs + 1
      if hsp.identity
        total_number_of_matches += hsp.identity
      elsif hsp.mismatch_count
        total_number_of_matches = hsp.align_len - hsp.mismatch_count
      end
    end
    return total_length_of_hsps, total_number_of_matches
  end
  # A method to retrieve a sequence from a blast datase using the fastacmd
  # @param [Hash] options The options for the method
  # @option options [String] :blast_database The path and name of the local blast database
  # @option options [String] :search_string The term to be searched for on the database, often the id or accession
  # @option options [String] :fastacmd_dir The directory in which to find the fastacmd executable
  def retrieve_sequence(options)
    default_options = {
      :fastacmd_dir => "",
      :protein_option =>"G" 
    }
    options.reverse_merge!(default_options)

    begin
      `#{options[:fastacmd_dir]}fastacmd 2>/dev/null`
    rescue
      puts "Can not find fastacmd at #{options[:fastacmd_dir]}fastacmd"
      raise
    end
    `#{options[:fastacmd_dir]}fastacmd -d #{options[:blast_database]} -p #{options[:protein_option]} -s #{options[:search_string]}`
  end
  # A method top retrieve the accession a hit from an ncbi blast search
  # @param [Bio::Blast::Report::Hit] blast_hit A blast hit object from which the accession will be
  # extracted
  def extract_accession_from_blast_hit(blast_hit)
    accession = blast_hit.definition.split("|")[3]
    accession.sub!(/\..+$/, "") # remove version number
  end
  # A method top retrieve a rich sequence object based on a hit from an ncbi blast search
  # @param [Bio::Blast::Report::Hit] blast_hit A blast hit object from which the accession will be
  # extracted and this used to retrieve the embl object
  def retrieve_ncbi_blast_hit_sequence(blast_hit)
    accession = extract_accession_from_blast_hit(blast_hit)
    retrieve_embl_object(accession)
  end
  # A method to retrieve an embl object by its accession
  # @param [String] accession The accession that will be used to retrieve the embl object
  def retrieve_embl_object(accession)
    server = Bio::Fetch.new('http://www.ebi.ac.uk/cgi-bin/dbfetch')
    embl_object = nil
    begin
      embl_text = server.fetch('embl', accession)
      embl_object = Bio::EMBL.new(embl_text)
    rescue
      embl_object = nil
    end
    return embl_object
  end
  # A method to retrieve an genbank refseq object by its accession
  # @param [String] accession The accession that will be used to retrieve the embl object
  def retrieve_refseq_object(accession)
    server = Bio::Fetch.new('http://www.ebi.ac.uk/cgi-bin/dbfetch')
    genbank_object = nil
    begin
      genbank_text = server.fetch('refseq', accession)
      genbank_object = Bio::GenBank.new(genbank_text)
    rescue
      genbank_object = nil
    end
    return genbank_object
  end

  # a method to return a local blast database friendly id from a rich sequence object
  # @param [Bio::(Sequence|Genbank::EMBL)] rich_sequence_object The rich sequence object from which the id will be generated
  def create_local_blast_database_id(rich_sequence_object)
    # database_id = rich_sequence_object.definition.sub(/#{rich_sequence_object.entry_id}\s+/,'').gsub(/\s/,'_').gsub(/,/,'')
    database_id = rich_sequence_object.entry_id
  end
  # a method to retrieve the sequence and qualifiers from the cds that matches a hit location
  # within a rich sequence object
  # @param [Bio::Blast::Report::Hit] hit The hit whose target start and end will be used to find the matching cds
  # @param [Bio::(Genbank|EMBL)] rich_sequence_object The rich sequence object from which to extract the cds details
  def retrieve_cds_hit_details_from_rich_sequence_object(hit,rich_sequence_object)
    if hit.target_start < hit.target_end
      earliest_location_in_hit = hit.target_start
    else
      earliest_location_in_hit = hit.target_end
    end
    rich_sequence_object.each_cds do |cds|
      location = cds.locations.first
      if earliest_location_in_hit >= location.from && earliest_location_in_hit < location.to
        cds_sequence = rich_sequence_object.to_biosequence.splice(cds.locations.to_s).to_s
        return cds_sequence,cds.qualifiers
      end
    end
    nil
  end
end
