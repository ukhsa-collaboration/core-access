##
# A set of methods to predict genes using glimmer3 and process the output
# 
# @author Anthony Underwood
module Glimmer
  require 'utils/hash_reverse_merge'
  require 'tempfile'
  require 'bioutils/rich_sequence_utils'

  require 'rubygems'
  require 'bio'

  # Predict genes using glimmer using the iterative approach that builds a training set from the long orfs of the input fasta
  # sequence and then predicts genes using the training set.
  # This will method will produce a predict file whose path is (glimmer_predict_filepath) will be returned
  # @param [Hash] options The options
  # @option options [String] :input_sequence_path The path to the input sequence (No spaces allowed)
  # @option options [Symbol] :input_format The format of the input sequence (can be guessed from input_sequence_path option)
  # @option options [String] :output_dir The output dir (by default this will be .)
  # @option option [String] :glimmer_dir_path the path to the glimmer directory (default is "/usr/local/glimmer/")
  # @option options [String] :rich_sequence_training_path The path to the file to be used for training. A rich sequence from which the cds can be extracted
  # @option options [String] cds_and_coords_training_path The path and prefix to the 3 files to be used for training e.g /home/user/sequences/GenomeA
  # will mean that in the directort/home/user/sequences there are 2 files GenomeA.fasta, GenomeA.train and GenomeA.coords where GenomeA.fasta
  # is the fasta sequence of the genome used for training, GenomeA.train is a multifasta of the coding sequences used for training and
  # GenomeA.coords is a coord file in the format gene_number\tstart\tend (ensuring direction for reverse complemented genes)
  # @option options [Boolean] :suppress_messages Suppress display of the verbose progress messages from glimmer or not (false by default)
  # @return [String] The path to the glimmer predict file produced
  def predict_genes_using_glimmer(options = {})
    default_options = {
      :output_dir => ".",
      :glimmer_dir_path => "/usr/local/glimmer",
      :rich_sequence_training_path => nil,
      :cds_and_coords_training_path => nil,
      :suppress_messages => false
    }
    options.reverse_merge!(default_options)
    options[:input_format] = guess_sequence_format(options[:input_sequence_path]) if options[:input_format].nil?

    case options[:input_format] # generate temp fasta file if input not fasta format
    when :genbank, :embl
      tmpfile_object = Tempfile.new('temp')
      temp_fasta_file = File.new(tmpfile_object.path,'w')
      case options[:input_format]
      when :genbank
        input_sequence = Bio::FlatFile.open(Bio::GenBank, options[:input_sequence_path]).next_entry
      when :embl
        input_sequence = Bio::FlatFile.open(Bio::EMBL, options[:input_sequence_path]).next_entry
      end
      temp_fasta_file.puts Bio::Sequence.new(input_sequence.seq).output(:fasta)
      temp_fasta_file.close
    end

    glimmer_predict_filename = File.basename(options[:rich_sequence_training_path], File.extname(options[:rich_sequence_training_path])) + "_glimmer"
    Dir.chdir(options[:output_dir]) # change to output dir in prepration for writing out glimmer files
    if options[:rich_sequence_training_path] || options[:cds_and_coords_training_path]
      if options[:rich_sequence_training_path] # make the necessary train and coords files derived from the cds sequences derived from a supplied rich sequence
        make_train_and_coord_files(:rich_sequence_training_path => options[:rich_sequence_training_path],
        :prefix => glimmer_predict_filename)
      end
      # set the path to the train and coord files depending on whether supplying a rich sequence file or just the coords and train files
      if options[:rich_sequence_training_path]
        cds_and_coords_training_path = "#{options[:output_dir]}/#{glimmer_predict_filename}"
      elsif options[:cds_and_coords_training_path]
        cds_and_coords_training_path = options[:cds_and_coords_training_path]
      end
      
      case options[:input_format]
      when :genbank, :embl
        predict_using_train_and_coords_file(:glimmer_dir_path => options[:glimmer_dir_path],
        :cds_and_coords_training_path => cds_and_coords_training_path,
        :input_sequence_path => tmpfile_object.path,
        :suppress_messages => options[:suppress_messages],
        :glimmer_predict_filename => glimmer_predict_filename)
      when :fasta
        predict_using_train_and_coords_file(:glimmer_dir_path => options[:glimmer_dir_path],
        :cds_and_coords_training_path => cds_and_coords_training_path,
        :input_sequence_path => options[:input_sequence_path],
        :suppress_messages => options[:suppress_messages],
        :glimmer_predict_filename => glimmer_predict_filename)
      else
        raise "Unrecognised sequence input format must be :fasta, :genbank or :embl"
      end
    else # predict based on the long orfs derived from the supplied fasta sequence
      case options[:input_format]
      when :genbank, :embl
        predict_using_iterated_glimmer(:glimmer_dir_path => options[:glimmer_dir_path],
        :input_sequence_path => tmpfile_object.path,
        :glimmer_predict_filename => glimmer_predict_filename,
        :suppress_messages => options[:suppress_messages])
      when :fasta
        predict_using_iterated_glimmer(:glimmer_dir_path => options[:glimmer_dir_path],
        :input_sequence_path => options[:input_sequence_path],
        :glimmer_predict_filename => glimmer_predict_filename,
        :suppress_messages => options[:suppress_messages])
      else
        raise "Unrecognised sequence input format must be :fasta, :genbak or :embl"
      end
    end
    glimmer_predict_filepath = "#{options[:output_dir]}/#{glimmer_predict_filename}.predict"
    return glimmer_predict_filepath
  end

  # A Method to predict genes using the iterated method where the sequence provided is fasta and the long orfs provide the training set
  # @param [hash] options The options
  # @option option [String] :glimmer_dir_path the path to the glimmer directory (default is "/usr/local/glimmer/")
  # @option options [String] :input_sequence_path The path to the fasta file that will be used as input
  # @option options [String] :glimmer_predict_filename The prefix for the glimmer predict file
  # @option options [Boolean] :suppress_messages Suppress display of the verbose progress messages from glimmer or not (false by default)
  def predict_using_iterated_glimmer(options)
    default_options = {
      :suppress_messages => false,
      :glimmer_dir_path => "/usr/local/glimmer",
    }
    options.reverse_merge!(default_options)
    raise ArgumentError, "the glimmer script can not be found at #{options[:glimmer_dir_path]}/scripts/" unless File.exists?("#{options[:glimmer_dir_path]}/scripts/g3-iterated.csh")
    command = "#{options[:glimmer_dir_path]}/scripts/g3-iterated.csh #{options[:input_sequence_path]} #{options[:glimmer_predict_filename]}"
    command += " 2>/dev/null" if options[:suppress_messages]
    `#{command}`
  end

  # # A Method to make the train and coords file when predicting genes via the non-iterated approach
  # @param [hash] options The options
  # @option options [String] :rich_sequence_training_path The path to the file to be used for training. A rich sequence from which the cds can be extracted
  # @option options [Symbol] :rich_sequence_format The format of the file to be used for training. Will be guessed if not suppplied
  # @option options [String] :prefix Prefix for output files
  def make_train_and_coord_files(options)
    options[:rich_sequence_format] = guess_sequence_format(options[:rich_sequence_training_path]) if options[:rich_sequence_format].nil?
    case options[:rich_sequence_format]
    when :genbank
      input_sequence_object = Bio::FlatFile.open(Bio::GenBank, options[:rich_sequence_training_path]).next_entry
    when :embl
      input_sequence_object = Bio::FlatFile.open(Bio::EMBL, options[:rich_sequence_training_path]).next_entry
    end
    input_sequence = input_sequence_object.seq
    
    #need fasta file during the glimmer training process
    fasta_ouput_file =  File.open("#{options[:prefix]}.fasta", "w")
    fasta_ouput_file.puts Bio::Sequence.new(input_sequence.seq).output(:fasta)
    fasta_ouput_file.close
    
    train_output_file = File.open("#{options[:prefix]}.train", "w")
    coords_output_file = File.open("#{options[:prefix]}.coords", "w")
    cds_count = 0
    input_sequence_object.each_cds do |cds|
      cds_count += 1
      train_output_file.puts "> gene #{cds_count}"
      coords_output_file.print "#{cds_count}\t"

      location =  cds.locations.first
      if location.strand == -1 #negative strand
        coords_output_file.puts "#{location.to}\t#{location.from}"
      else
        coords_output_file.puts "#{location.from}\t#{location.to}"
      end

      location_string = location.range.to_s
      if location.strand == -1 #negative strand
        location_string = "complement(#{location_string})"
      end

      cds_sequence = input_sequence.splice(location_string)

      train_output_file.puts cds_sequence
    end
    train_output_file.close
    coords_output_file.close
  end

  # A Method to predict genes based on a training set derived from cds and supplied as train and coords files. N.B Need to make sure formatdb is in the PATH
  # @param [hash] options The options
  # @option option [String] :glimmer_dir_path the path to the glimmer directory (default is "/usr/local/glimmer/")
  # @option options [String] cds_and_coords_training_path The path and prefix to the 3 files to be used for training e.g /home/user/sequences/GenomeA
  # will mean that in the directort/home/user/sequences there are 2 files GenomeA.fasta, GenomeA.train and GenomeA.coords where GenomeA.fasta
  # is the fasta sequence of the genome used for training, GenomeA.train is a multifasta of the coding sequences used for training and
  # GenomeA.coords is a coord file in the format gene_number\tstart\tend (ensuring direction for reverse complemented genes)
  # @option options [Boolean] :suppress_messages Suppress display of the verbose progress messages from glimmer or not (false by default)
  # @option options [String] :input_sequence_path The path to the fasta file that will be used as input
  # @option options [String] :glimmer_predict_filename The prefix for the glimmer predict filename
  def predict_using_train_and_coords_file(options)
    prefix = File.basename(options[:cds_and_coords_training_path])
    startuse = create_model(:glimmer_dir_path => options[:glimmer_dir_path], :cds_and_coords_training_path => options[:cds_and_coords_training_path], :prefix => prefix,  :suppress_messages => options[:suppress_messages])
    run_glimmer_using_model(:glimmer_dir_path => options[:glimmer_dir_path], :input_sequence_path => options[:input_sequence_path], :glimmer_predict_filename => options[:glimmer_predict_filename], :startuse => startuse, :prefix => prefix, :suppress_messages => options[:suppress_messages])
  end
  # A Method to predict a model for glimmer
  # @param [hash] options The options
  # @option option [String] :glimmer_dir_path the path to the glimmer directory (default is "/usr/local/glimmer/")
  # @option options [String] cds_and_coords_training_path The path and prefix to the 3 files to be used for training e.g /home/user/sequences/GenomeA
  # will mean that in the directort/home/user/sequences there are 2 files GenomeA.fasta, GenomeA.train and GenomeA.coords where GenomeA.fasta
  # is the fasta sequence of the genome used for training, GenomeA.train is a multifasta of the coding sequences used for training and
  # GenomeA.coords is a coord file in the format gene_number\tstart\tend (ensuring direction for reverse complemented genes)
  # @option option [String] prefix prefix of the files to be created
  def create_model(options)
    default_options = {
      :glimmer_dir_path => "/usr/local/glimmer",
    }
    options.reverse_merge!(default_options)
    # build a model
    `#{options[:glimmer_dir_path]}/bin/build-icm -r #{options[:prefix]}.icm < #{options[:cds_and_coords_training_path]}.train`
    # make  upstream info
    `#{options[:glimmer_dir_path]}/scripts/upstream-coords.awk 25 0 #{options[:cds_and_coords_training_path]}.coords | #{options[:glimmer_dir_path]}/bin/extract #{options[:cds_and_coords_training_path]}.fasta - > #{options[:prefix]}.upstream`
    # Make a PWM of the upstream regions
    command = "elph #{options[:prefix]}.upstream LEN=6"
    command += " 2>/dev/null " if options[:suppress_messages]
    command += "| #{options[:glimmer_dir_path]}/scripts/get-motif-counts.awk > #{options[:prefix]}.motif"
    `#{command}`
    # set the startuse variable
    startuse = `#{options[:glimmer_dir_path]}/bin/start-codon-distrib -3 #{options[:cds_and_coords_training_path]}.fasta #{options[:cds_and_coords_training_path]}.coords`.chomp
    return startuse
  end

  # # Run gliimer3 using training model
  # @param [hash] options The options
  # @option option [String] :glimmer_dir_path the path to the glimmer directory (default is "/usr/local/glimmer/")
  # @option option [String] :startuse A string of the start codon usages
  # @option options [Boolean] :suppress_messages Suppress display of the verbose progress messages from glimmer or not (false by default)
  # @option options [String] :input_sequence_path The path to the fasta file that will be used as input
  # @option options [String] :glimmer_predict_filename The prefix for the glimmer predict filename
  def run_glimmer_using_model(options)
    default_options = {
      :suppress_messages => false,
      :glimmer_dir_path => "/usr/local/glimmer",
    }
    options.reverse_merge!(default_options)
    unless options[:glimmer_predict_filename]
      options[:glimmer_predict_filename] = File.basename(options[:input_sequence_path], File.extname(options[:input_sequence_path])) + "_glimmer"
    end
    unless options[:startuse]
      options[:startuse] = `#{options[:glimmer_dir_path]}/bin/start-codon-distrib -3 #{options[:prefix]}.fasta #{options[:prefix]}.coords`.chomp
    end
    # run glimmer to predict genes
    command = "#{options[:glimmer_dir_path]}/bin/glimmer3 -o50 -g150 -t30 -l -b #{options[:prefix]}.motif -P #{options[:startuse]} #{options[:input_sequence_path]} #{options[:prefix]}.icm #{options[:glimmer_predict_filename]}"
    # command = "/usr/local/glimmer-mg/bin/glimmer-mg -i -o 50 -g 150 -b #{prefix}.motif -m {prefix}.icm #{options[:input_sequence_path]}  #{options[:glimmer_predict_filename]}"
    command += " 2>/dev/null " if options[:suppress_messages]
    `#{command}`
  end


  # This method will generate a genbank file that contains CDS features based on glimmer gene predictions. These will have no
  # annotation other than that they are coding sequences i.e no further functional prediction.
  # @param [Hash] options The options used in teh method
  # @option options [String] :glimmer_predict_file The path to the glimmer predict file
  # @option options [String] :input_sequence_path The path to the sequence file to be annotated
  # @option options [Symbol] :input_format The input format of the sequence to be annotated
  # @option options [String] :output_dir path to the output directory
  # @option options [Symbol] :output_format The format in which the rich sequence will be output
  # @return [String] The path to the output file
  def glimmer_prediction_to_rich_sequence_file(options = {})
    default_options = {
      :output_dir => ".",
      :output_format => :genbank
    }
    options.reverse_merge!(default_options)

    predicted_genes = File.open(options[:glimmer_predict_file])
     options[:input_format] = guess_sequence_format(options[:input_sequence_path]) if options[:input_format].nil?
    case options[:input_format]
    when :fasta
      input_sequences = Bio::FlatFile.open(Bio::FastaFormat, options[:input_sequence_path])
    when :genbank
      input_sequences = Bio::FlatFile.open(Bio::GenBank, options[:input_sequence_path])
    when :embl
      input_sequences = Bio::FlatFile.open(Bio::EMBL, options[:input_sequence_path])
    end

    bio_sequence = nil
    features = nil
    input_sequence = nil
   
    #get ready to write out to a file
    output_filename = File.basename(options[:input_sequence_path], File.extname(options[:input_sequence_path])) + "_glimmer_genes"
    case options[:output_format]
    when :genbank
      output_filename += ".gbk"
    when :embl
      output_filename += ".embl"
    end
    File.delete("#{options[:output_dir]}/#{output_filename}") if File.exists?("#{options[:output_dir]}/#{output_filename}")# clear old data
    output_file = File.open("#{options[:output_dir]}/#{output_filename}", "a")
  
    predicted_genes.each do |line|
      if line =~ /^>/#header line
        unless input_sequence.nil?
          case options[:input_format]
          when :genbank, :embl
            bio_sequence.features = features + input_sequence.features # incorporate existing features
          when :fasta
            bio_sequence.features = features
          end
          case options[:input_format]
          when :genbank, :embl
            bio_sequence.features = features + input_sequence.features # incorporate existing features
          when :fasta
            bio_sequence.features = features
          end
          output_file.puts bio_sequence.output(options[:output_format])
        end
        input_sequence = input_sequences.next_entry
        bio_sequence = Bio::Sequence.new(input_sequence.seq)
        bio_sequence.entry_id = input_sequence.entry_id
        bio_sequence.na
        bio_sequence.definition = input_sequence.definition
        bio_sequence.date_created = "#{Date.today.day}-#{Date::ABBR_MONTHNAMES[Date.today.month].upcase}-#{Date.today.year}"
        features = []
      else
        orf_name, start_pos , end_pos, frame, score = line.split(/\s+/)
        if frame =~ /\+/ # check for frame and swap start and end if necessary
          feature = Bio::Feature.new('CDS',"#{start_pos}..#{end_pos}")
        else
          feature = Bio::Feature.new('CDS',"complement(#{end_pos}..#{start_pos})")
        end
        features << feature
      end
    end
  
    case options[:input_format]
    when :genbank, :embl
      bio_sequence.features = features + input_sequence.features # incorporate existing features
    when :fasta
      bio_sequence.features = features
    end
    output_file.puts bio_sequence.output(options[:output_format])
    
    output_file.close
    return "#{options[:output_dir]}/#{output_filename}"
  end
end
