class GenomeReciprocalHitAnnotator
  require 'utils/method_argument_parser'
  require 'utils/hash_reverse_merge'
  
  require 'bioutils/glimmer'
  require 'bioutils/blast'
  require 'bioutils/reciprocal_hit'
  extend Glimmer
  extend Blast
  extend ReciprocalHit
  class << self
    # a method to predict genes in a fasta sequence and annotate from reference genome(s) at first and then versus ncbi if no hits to the reference
    # @param [Hash] options The options for the method
    # @option options [String] :sequence_to_be_annotated_path The path to the sequence that will be annotated (in fasta format)
    # @option options [String] :training_sequence_path The path to the rich sequence that wil be used as a training set for
    # gene prediction. If nil then the long orfs from the sequence to be annotated will be used
    # @option options [Array] :reference_genomes An array of path strings each one pointing to a reference genome that will be used for annotation
    # @option options [String] :output_dir The path to location of a directory where the output files will be written
    # @option options [String] :path_to_blast_executable The path to the blast executable, default is "/usr/local/blast/bin/blastall",
    # @option options [String] :fastacmd_dir The directory in which to find the fastacmd executable (default is "/usr/local/blast/bin/")
    # @option options [String] :e_cutoff The expect cutoff to use when doing blast searches (a smaller number here e.g 1e-20 this can speed
    # up the blast searches)
    # @option options [Boolean] :accept_reciprocal_hits_with_multiple_hits_incl_query This determines whether to accept remote blast hits where
    # the reciprocal hit does has multiple hits (including the original query) (i.e possible gene family). Default is true
    # @option options [Boolean] :suppress_messages Whether to suppress messages when generating glimmer files etc (default is true)
    # @option options [Symbol] :output_format The type of output to be generated (:genbank or :embl). The default output is genbank
    # These should be in genbank or embl formats
    # @option options [Boolean] :annotate_by_remote_blast Whether to annotate by blasting against the NCBI nr database if a local hit isn't found. Default is true 
    # @option options [Boolean] :annotate_vs_local_microbial_genomes Whether to annotate by blasting against the local database of all viral & bacterial genomes if no hit to reference genomes. Default is true.
    # @option options [String] :reference_blast_program Whether to use blastn or blastp for reference genome blast. Default is blastn
    # @option options [String] :reference_percent_identity_cutoff Percent id cut off for significant hit when doing local blast versus ref genomes
    # @option options [String] :reference_minimum_hit_length minimum hit length when doing local blast versus ref genomes
    # @option options [String] :microbial_genomes_blast_program Whether to use blastn or blastp for microbial genomes blast. Default is blastn
    # @option options [String] :microbial_genomes_percent_identity_cutoff Percent id cut off for significant hit when doing local blast versus microbial genomes
    # @option options [String] :microbial_genomes_minimum_hit_length minimum hit length when doing local blast versus microbial genomes
    # @option options [String] :ncbi_blast_program Whether to use blastn or blastp for NCBI blast. Default is blastn
    # @option options [String] :ncbi_percent_identity_cutoff Percent id cut off for significant hit when doing remote blast versus NCBI nr
    # @option options [String] :ncbi_minimum_hit_length minimum hit length when doing remote blast versus NCBI nr
    # @options options [Boolean] :display_metrics Whether or not to show results of significant hit testing
    # @options options [Boolean] :display_query_sequence Whether or not to display the query sequence. Useful for debugging. Default false
    # @options options [Boolean] :debug Whether to write a log file. Default false
    # @options options [String] :log_file_location Location of the log file
    # @options options [Integer] :parallel_processors If not nil(default) the number of processors to use in parrallel annotation
    def generate_annotation_by_reciprocal_hits(options)
      default_options = {
        :path_to_blast_executable => "/usr/local/blast/bin/blastall",
        :fastacmd_dir => "/usr/local/blast/bin/",
        :e_cutoff => "1e-20",
        :accept_reciprocal_hits_with_multiple_hits_incl_query => true,
        :suppress_messages => true,
        :output_format => :genbank,
        :parallel_processors => nil,
        :annotate_by_remote_blast => true,
        :annotate_vs_local_microbial_genomes => true,
        :local_microbial_blast_DB => "/Volumes/DataRAID/blast_databases/microbial_genomes",
        :reference_blast_program => "blastn",
        :reference_percent_identity_cutoff => 90,
        :reference_minimum_hit_length => 85,
        :microbial_genomes_blast_program => "blastn",
        :microbial_genomes_percent_identity_cutoff => 90,
        :microbial_genomes_minimum_hit_length => 85,
        :ncbi_blast_program => "blastn",
        :ncbi_percent_identity_cutoff => 90,
        :ncbi_minimum_hit_length => 85,
        :display_metrics => false,
        :display_query_sequence => false,
        :debug => false,
        :log_file_location => Dir.pwd
      }
      options.reverse_merge!(default_options)
      
      options = MethodArgumentParser::Parser.check_options options  do
        option :sequence_to_be_annotated_path, :required => true, :type => :string
        # option :reference_genomes,  :required => true, :type => :array
      end

      Dir.chdir(options[:output_dir]) if options[:output_dir]

      if options[:debug]
        debug = File.open("#{options[:log_file_location]}/annotation_debug.log", "a")
      end
      
      # predict genes using glimmer and obtain path to predict file      
      if File.exists?(File.basename(options[:sequence_to_be_annotated_path], File.extname(options[:sequence_to_be_annotated_path])) + "_glimmer.predict")
        predict_file = File.basename(options[:sequence_to_be_annotated_path], File.extname(options[:sequence_to_be_annotated_path])) + "_glimmer.predict"
      else
        puts "Predicting genes using glimmer ......" 
        predict_file = predict_genes_using_glimmer(:input_sequence_path  => options[:sequence_to_be_annotated_path],
                                                 :rich_sequence_training_path => options[:training_sequence_path],
                                                 :suppress_messages => false)
      end
      # convert the glimmer prediction to a gebank file with the gene CDS marked
      puts "Converting glimmer prediction to a genbank file ....."
      glimmer_genbank_file = glimmer_prediction_to_rich_sequence_file(:glimmer_predict_file => predict_file,
                                                                       :input_sequence_path  => options[:sequence_to_be_annotated_path])


      # now start the annotation
      puts "Starting annotation....."
      puts "Making Blast databases ....."
      # databases from glimmer predicted genes for reverse blast in the reciprocal blast
      if options[:reference_blast_program] == "blastp" || options[:microbial_genomes_blast_program] == "blastp" || options[:ncbi_blast_program] == "blastp"
        query_sequence_database_name = blast_database_from_rich_sequences(:input_sequences => [glimmer_genbank_file], :protein => true)
      end
      if options[:reference_blast_program] == "blastn" || options[:microbial_genomes_blast_program] == "blastn" || options[:ncbi_blast_program] == "blastn"
        query_sequence_database_name = blast_database_from_rich_sequences(:input_sequences => [glimmer_genbank_file], :protein => false)
      end
      # databases from sequences that contain annotation
      if options[:reference_database_name]
        reference_sequences_database_name = options[:reference_database_name]
      else
        options[:reference_genomes].size > 1 ? ref_genome_database_name = "Reference_Genomes" : ref_genome_database_name = nil
        
        if options[:reference_blast_program] == "blastp" || options[:microbial_genomes_blast_program] == "blastp" || options[:ncbi_blast_program] == "blastp"
          reference_sequences_database_name = blast_database_from_rich_sequences(:input_sequences => options[:reference_genomes], :database_labels => :full_annotation, :database_name => ref_genome_database_name, :protein => true)
        end
        if options[:reference_blast_program] == "blastn" || options[:microbial_genomes_blast_program] == "blastn" || options[:ncbi_blast_program] == "blastn"
          reference_sequences_database_name = blast_database_from_rich_sequences(:input_sequences => options[:reference_genomes], :database_labels => :full_annotation, :database_name => ref_genome_database_name, :protein => false)
        end
        
      end

      #get output file ready
      output_filename = File.dirname(options[:sequence_to_be_annotated_path]) + "/" + File.basename(options[:sequence_to_be_annotated_path], File.extname(options[:sequence_to_be_annotated_path])) + "-annotated.gbk"
      File.delete(output_filename) if File.exists?(output_filename)# clear old data
      output_file = File.open(output_filename, "w")
      
      input_sequence_with_precicted_genes = Bio::FlatFile.open(Bio::GenBank, glimmer_genbank_file)
      if options[:parallel_processors]
        require 'forkoff'
        
        combined_annotations = Hash.new # {sequence_def1 => [cds, cds, cds ....], sequence_def2 => [cds, cds, cds ...]}
        input_sequence_objects = Hash.new # {sequence_def1 => sequence_object1, sequence_def2 => sequence_object2, ...}
        input_biosequences  = Hash.new # {sequence_def1 => biosequence1, sequence_def2 => biosequence2, ...}
        combined_coding_sequences = Array.new # a big array of arrays [[sequence_def1, 1, cds], [sequence_def1, 2, cds], ....]
        # loop through all contigs
        while input_sequence_with_precicted_genes_object = input_sequence_with_precicted_genes.next_entry
        # collect all the CDS features into Array
          coding_sequences = input_sequence_with_precicted_genes_object.features.select{|f| f.feature == 'CDS'}
          if coding_sequences.size > 0
            cds_number = 0
            combined_coding_sequences += coding_sequences.map{|coding_sequence| [input_sequence_with_precicted_genes_object.definition, cds_number += 1, coding_sequence]}
          end
          # record contig in hash for later use
          input_sequence_objects[input_sequence_with_precicted_genes_object.definition] = input_sequence_with_precicted_genes_object
          input_biosequences[input_sequence_with_precicted_genes_object.definition] = input_sequence_with_precicted_genes_object.seq
          combined_annotations[input_sequence_with_precicted_genes_object.definition] = Array.new
        end
            
        
        # work out slice size, put into chunks of 10 so a fork can only slow down the overall process so much
        if (combined_coding_sequences.size.to_f/options[:parallel_processors]).ceil < 10
          slice_size = (combined_coding_sequences.size.to_f/options[:parallel_processors]).ceil
        else
          slice_size = 10
        end
        
        # make an array of array of each slice
        combined_coding_sequences_slices = Array.new
        cds_counter = 0
        combined_coding_sequences.each_slice(slice_size) do |combined_coding_sequences_slice|
          combined_coding_sequences_slices << combined_coding_sequences_slice
        end
        
        # parallel loop starts here
        slice_number = 0
        collected_features = combined_coding_sequences_slices.forkoff :processes => options[:parallel_processors], :strategy => :file do |*combined_coding_sequences_slice|
          slice_annotated_features = Array.new
          combined_coding_sequences_slice.each do |cds_details|
            cds = cds_details[2]
            cds_number = cds_details[1]
            sequence_definition = cds_details[0]
            cds_location_string = cds.locations.to_s
            cds_sequence = input_biosequences[sequence_definition].splice(cds_location_string)
            # annotated_cds = annotate_feature(cds, cds_number, input_biosequences[sequence_definition], input_sequence_objects[sequence_definition], query_sequence_database_name, reference_sequences_database_name, options, debug)
            if options[:debug]
              puts "#{File.basename(Dir.pwd)} #{sequence_definition}: #{cds_number} (#{cds_location_string})"
              debug.puts "#{File.basename(Dir.pwd)} #{sequence_definition}: #{cds_number} (#{cds_location_string})"
              debug.flush
            end
            annotated_cds = annotate_feature(cds, options.merge(:entry_id => "#{input_sequence_objects[sequence_definition].entry_id}-#{cds_number}", :sequence => cds_sequence, :query_sequence_database_name => query_sequence_database_name, :reference_sequences_database_name => reference_sequences_database_name, :debug_file => debug))
            slice_annotated_features << [sequence_definition, annotated_cds]
          end
          slice_annotated_features
        end
        collected_features.each do |slice_annotated_features|
          slice_annotated_features.each do |sequence_definition, annotated_feature|
            combined_annotations[sequence_definition] << annotated_feature
          end
        end
        input_sequence_objects.each_key do |sequence_definition|
          input_sequence_with_precicted_genes_object = input_sequence_objects[sequence_definition]
          input_biosequence = input_biosequences[sequence_definition]

          # get sequence object ready for annotated sequence
          annotated_sequence = Bio::Sequence.new(input_biosequence)
          annotated_sequence.entry_id = input_sequence_with_precicted_genes_object.entry_id
          annotated_sequence.na
          annotated_sequence.definition = input_sequence_with_precicted_genes_object.definition
          annotated_sequence.date_created = "#{Date.today.day}-#{Date::ABBR_MONTHNAMES[Date.today.month].upcase}-#{Date.today.year}"
          annotated_sequence.molecule_type = "DNA"
          annotated_sequence.topology = "linear"
          annotated_sequence.division = "PRO"
          annotated_sequence.features = combined_annotations[input_sequence_with_precicted_genes_object.definition]

          #write out file
          output_file.puts annotated_sequence.output(options[:output_format])
        end
      else # non parallel
        while input_sequence_with_precicted_genes_object = input_sequence_with_precicted_genes.next_entry
          puts "Annotating #{input_sequence_with_precicted_genes_object.definition}"
          input_biosequence = input_sequence_with_precicted_genes_object.seq

          # get sequence object ready for annotated sequence
          annotated_sequence = Bio::Sequence.new(input_biosequence)
          annotated_sequence.entry_id = input_sequence_with_precicted_genes_object.entry_id
          annotated_sequence.na
          annotated_sequence.definition = input_sequence_with_precicted_genes_object.definition
          annotated_sequence.date_created = "#{Date.today.day}-#{Date::ABBR_MONTHNAMES[Date.today.month].upcase}-#{Date.today.year}"
          annotated_sequence.molecule_type = "DNA"
          annotated_sequence.topology = "linear"
          annotated_sequence.division = "PRO"
          annotated_features = input_sequence_with_precicted_genes_object.features.dup
      
          cds_counter = 0
          embl_objects = {}

          input_sequence_with_precicted_genes_object.each_cds do |cds|
            cds_counter += 1
            # next  if cds_counter == 1 # debugging line
            # break if cds_counter > 50 # debugging line
            # print info for counter

            # annotate_feature(cds, cds_counter, input_biosequence, input_sequence_with_precicted_genes_object, query_sequence_database_name, reference_sequences_database_name, options, debug, embl_objects)
            annotate_feature(cds,options.merge(:debug_file => debug))
          end
          annotated_sequence.features = annotated_features

          #write out file
          output_file.puts annotated_sequence.output(options[:output_format])
        end
      end
      output_file.close
      if options[:debug]
        debug.close
      end
    end
    
    # method to add an array of qualifiers to a feature
    # @param [Bio::Feature] feature The feature which will have the qualifiers appended to it
    # @param [Array] qualifiers An array of two element arrays where the first element in the array is the qualifier type
    # e.g 'note' and the second element is the value e.g 'possible virulence factor'
    def add_qualifiers_to_feature(feature, qualifiers)
      qualifiers.each do |qualifier|
        qualifier_type = qualifier[0]
        qualifier_value = qualifier[1]
        feature.append(Bio::Feature::Qualifier.new(qualifier_type, qualifier_value)) unless qualifier_type =~ /(translation|transl_table)/i # append annotation apart from translation
      end
    end

    def annotate_feature(feature,options)
      options = MethodArgumentParser::Parser.check_options options  do
        option :sequence, :required => true
        option :entry_id, :required => true
      end

      biosequence = Bio::Sequence.new(options[:sequence])
      biosequence.entry_id = options[:entry_id]

      reciprocal_hit_qualifiers = annotate_sequence(options.merge(:biosequence => biosequence))
      add_qualifiers_to_feature(feature, reciprocal_hit_qualifiers) unless reciprocal_hit_qualifiers.nil?
      return feature
    end
    
    def annotate_sequence(options)
      default_options = {
        :path_to_blast_executable => "/usr/local/blast/bin/blastall",
        :fastacmd_dir => "/usr/local/blast/bin/",
        :e_cutoff => "1e-20",
        :annotate_by_remote_blast => true,
        :annotate_vs_local_microbial_genomes => true,
        :embl_objects => nil
      }
      options.reverse_merge!(default_options)
      
      options = MethodArgumentParser::Parser.check_options options  do
        option :biosequence, :required => true
        option :query_sequence_database_name, :required => true, :type => :string
        option :reference_sequences_database_name, :required => true, :type => :string
        option :reference_blast_program, :required => true, :type => :string
        option :reference_percent_identity_cutoff, :required => true, :type => :integer
        option :reference_minimum_hit_length, :required => true, :type => :integer
        option :local_microbial_blast_DB, :required => true, :type => :string
      end

      nucleotide_query_sequence = options[:biosequence]

      protein_query_sequence = Bio::Sequence.new(nucleotide_query_sequence.translate(1,11))
      protein_query_sequence.entry_id = nucleotide_query_sequence.entry_id
      protein_query_sequence.aa 
      

      # try local blast first
      if options[:reference_blast_program] == "blastn"
        query_sequence = nucleotide_query_sequence
      else
        query_sequence = protein_query_sequence
      end
      if options[:display_query_sequence]
        puts query_sequence
        if options[:debug]
          options[:debug_file].puts query_sequence
          options[:debug_file].flush
        end
      end
      reciprocal_hit_found, reciprocal_hit_details, reciprocal_hit_qualifiers = find_local_reciprocal_hit(:query_sequence => query_sequence,:query_sequence_blast_database => options[:query_sequence_database_name], :subject_blast_database => options[:reference_sequences_database_name], :path_to_blast_executable => options[:path_to_blast_executable], :blast_program => options[:reference_blast_program], :blast_options => "-e #{options[:e_cutoff]} -F F -b 10 -v 10",:fastacmd_dir => options[:fastacmd_dir], :percent_identity_cutoff => options[:reference_percent_identity_cutoff], :minimum_hit_length => options[:reference_minimum_hit_length], :display_metrics => options[:display_metrics])
      if reciprocal_hit_found || (!reciprocal_hit_qualifiers.nil? && options[:accept_first_reciprocal_hit_containing_query])
        puts "Found hit to #{reciprocal_hit_details} by local blast"
        reciprocal_hit_qualifiers.unshift(["note", "Annotation derived by reciprocal #{options[:reference_blast_program]} analysis against reference sequences (percent_identity_cutoff: #{options[:reference_percent_identity_cutoff]}, minimum_hit_length: #{options[:reference_minimum_hit_length]}). Annotation from #{extract_organism_qualifier(reciprocal_hit_qualifiers)}"])
        return reciprocal_hit_qualifiers
      end
      # local blast failed try against all microbial genomes
      if options[:annotate_vs_local_microbial_genomes] && !reciprocal_hit_found
        # blast vs reference genomes failed, try vs local all-microbe database
        if reciprocal_hit_details !~ /^No/
          puts "local blast failed vs ref genomes: Matched #{reciprocal_hit_details} but found multiple hits with the reciprocal blast"
        else
          puts "local blast failed vs ref genomes: #{reciprocal_hit_details}"
        end
        if options[:microbial_genomes_percent_identity_cutoff]
          percent_identity_cutoff = options[:microbial_genomes_percent_identity_cutoff]
        else
          if options[:blast_program] == "blastp"
            percent_identity_cutoff = 70
          else
            percent_identity_cutoff = 90
          end
        end
        if options[:microbial_genomes_blast_program] == "blastn"
          query_sequence = nucleotide_query_sequence
        else
          query_sequence = protein_query_sequence
        end
        if options[:display_query_sequence]
          puts query_sequence
          if options[:debug]
            options[:debug_file].puts query_sequence
            options[:debug_file].flush
          end
        end
        reciprocal_hit_found, reciprocal_hit_details, reciprocal_hit_qualifiers = find_local_reciprocal_hit(:query_sequence => query_sequence,:query_sequence_blast_database => options[:query_sequence_database_name],:subject_blast_database => options[:local_microbial_blast_DB], :path_to_blast_executable => options[:path_to_blast_executable], :blast_program => options[:microbial_genomes_blast_program], :blast_options => "-e #{options[:e_cutoff]} -F F -b 10 -v 10",:percent_identity_cutoff => percent_identity_cutoff, :minimum_hit_length => options[:microbial_genomes_minimum_hit_length], :fastacmd_dir => options[:fastacmd_dir], :display_metrics => options[:display_metrics])
        if reciprocal_hit_found || (!reciprocal_hit_qualifiers.nil? && options[:accept_first_reciprocal_hit_containing_query])
          puts "Found hit to #{reciprocal_hit_details} by local blast against all microbial genomes"
          reciprocal_hit_qualifiers.unshift(["note", "Annotation derived by reciprocal #{options[:microbial_genomes_blast_program]} analysis against all microbial genome sequences (percent_identity_cutoff: #{percent_identity_cutoff}, minimum_hit_length: #{options[:microbial_genomes_minimum_hit_length]}). Annotation from #{extract_organism_qualifier(reciprocal_hit_qualifiers)}"])
          return reciprocal_hit_qualifiers
        elsif !reciprocal_hit_qualifiers.nil? && !options[:annotate_by_remote_blast] && options[:accept_reciprocal_hits_with_multiple_hits_incl_query]
          puts "Found hit to #{reciprocal_hit_details} by local blast against all microbial genomes"
          reciprocal_hit_qualifiers.unshift(["note", "N.B reciprocal blast analysis found that the reciprocal hit matched multiple sequences inclusing the query (possible gene family). Annotation derived by reciprocal #{options[:microbial_genomes_blast_program]} analysis against all microbial genome sequences (percent_identity_cutoff: #{percent_identity_cutoff}, minimum_hit_length: #{options[:microbial_genomes_minimum_hit_length]}). Annotation from #{extract_organism_qualifier(reciprocal_hit_qualifiers)}"])
          return reciprocal_hit_qualifiers
        end
      end
      # local blast vs ref genomes and microbial genomes failed
      if options[:annotate_by_remote_blast] && !reciprocal_hit_found
        if reciprocal_hit_details !~ /^No/
          puts "local blast failed vs microbial genomes: Matched #{reciprocal_hit_details} but found multiple hits with the reciprocal blast"
        else
          puts "local blast failed vs microbial genomes: #{reciprocal_hit_details}"
        end
        if options[:ncbi_percent_identity_cutoff]
          percent_identity_cutoff = options[:ncbi_percent_identity_cutoff]
        else
          if options[:blast_program] == "blastp"
            percent_identity_cutoff = 70
          else
            percent_identity_cutoff = 90
          end
        end
        if options[:ncbi_blast_program] == "blastn"
          query_sequence = nucleotide_query_sequence
        else
          query_sequence = protein_query_sequence
        end
        if options[:display_query_sequence]
          if options[:debug]
            options[:debug_file].puts query_sequence
            options[:debug_file].flush
          end
        end
        reciprocal_hit_found, reciprocal_hit_details, reciprocal_hit_qualifiers, embl_object = find_remote_reciprocal_hit(:query_sequence => query_sequence,:query_sequence_blast_database => options[:query_sequence_database_name],:blast_program => options[:ncbi_blast_program], :blast_options => {:expect => options[:e_cutoff]}, :percent_identity_cutoff => percent_identity_cutoff, :minimum_hit_length => options[:ncbi_minimum_hit_length],:embl_objects => options[:embl_objects], :display_metrics => options[:display_metrics])
        if reciprocal_hit_found 
          options[:embl_objects][embl_object.accession] = embl_object unless embl_object.nil? || options[:embl_objects].nil? || options[:embl_objects].has_key?(embl_object.accession)
          puts "Found hit to #{reciprocal_hit_details} by remote blast"
          reciprocal_hit_qualifiers.unshift(["note", "Annotation derived by reciprocal #{options[:ncbi_blast_program]} analysis against NCBI nr (percent_identity_cutoff: #{percent_identity_cutoff}, minimum_hit_length: #{options[:ncbi_minimum_hit_length]}). Annotation from #{extract_organism_qualifier(reciprocal_hit_qualifiers)}"])
          return reciprocal_hit_qualifiers
        elsif !reciprocal_hit_qualifiers.nil? && options[:accept_reciprocal_hits_with_multiple_hits_incl_query]
          puts "Found hit to #{reciprocal_hit_details} by remote blast but the reciprocal blast returned a sequence other than the query (possible gene family)"
          reciprocal_hit_qualifiers.unshift(["note", "N.B reciprocal blast analysis found that the reciprocal hit matched multiple sequences inclusing the query (possible gene family). Annotation derived by reciprocal #{options[:ncbi_blast_program]} analysis against NCBI nr (percent_identity_cutoff: #{percent_identity_cutoff}, minimum_hit_length: #{options[:ncbi_minimum_hit_length]}). Annotation from #{extract_organism_qualifier(reciprocal_hit_qualifiers)}"])
          return reciprocal_hit_qualifiers
        else
          puts reciprocal_hit_details
          return nil
        end
      end
    end
  end
end
