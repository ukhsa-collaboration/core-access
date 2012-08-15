module ClusterCreate
  gem "bio", "~> 1.4.2"
  require 'bio'
  require 'utils/hash_reverse_merge'
  require 'utils/method_argument_parser'

  # Create a multi fasta file containing the genes of all the input sequences. If an input sequence is in fasta format then the genes will be predicted using Glimmer. Glimmer can be trained with a training sequence or with an existing glimmer model or just using itereated glimmer. This depends on the options supplied (see below)
  # @param [Hash] options The options for the method
  # @option options [String] :root_folder the path and name of the folder where the results files will be written
  # @option options [String] :cds_multi_fasta_file The name of the multi fasta that will be created
  # @option options [Array] :sequence_files An array of filepaths listing each of the input sequence files
  # @option options [String] :training_sequence_path The path to a training sequence to be used in Glimmer  
  # @option options [String] :training_model_prefix The path to the training model if a training model has already been created
  # @option options [Boolean] :verbose Print verbose output

  def create_cds_multi_fasta_file(options)
    require 'bioutils/rich_sequence_utils'
    require 'bioutils/glimmer'
    extend Glimmer

    default_options = {
      :cds_multi_fasta_file => "cds_proteins.fas",
      :verbose => false
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :root_folder, :required => true, :type => :string
      option :cds_multi_fasta_file, :required => true, :type => :string
      option :sequence_files, :required => true, :type => :array

    end

    Dir.chdir(options[:root_folder])

    files_with_cds = Array.new # a list of files containing
    options[:sequence_files].each do |sequence_file|
      sequence_format = guess_sequence_format(sequence_file)
      if sequence_format == :fasta
        if options[:training_model_prefix]
          puts "Predicting genes for file #{sequence_file} using training model ...."
          run_glimmer_using_model(:input_sequence_path  => sequence_file, :prefix => options[:training_model_prefix],:glimmer_dir_path => options[:glimmer_dir], :suppress_messages => true)
          predict_file = File.basename(sequence_file, File.extname(sequence_file)) + "_glimmer.predict"
        elsif options[:training_sequence_path]
          model_file_prefix = File.basename(options[:training_sequence_path], File.extname(options[:training_sequence_path])) + "_glimmer"
          if File.exists?(model_file_prefix + ".icm")
            if options[:verbose]
              puts "Predicting genes for file #{sequence_file} using training model ...."
            else
              print "."
            end
            run_glimmer_using_model(:input_sequence_path  => sequence_file, :prefix => model_file_prefix,:glimmer_dir_path => options[:glimmer_dir], :suppress_messages => true)
            predict_file = File.basename(sequence_file, File.extname(sequence_file)) + "_glimmer.predict"
          else
            if options[:verbose]
              puts "Predicting genes for file #{sequence_file} using training sequence ...."
            else
              print "."
            end
            predict_file = predict_genes_using_glimmer(:input_sequence_path  => sequence_file,
                                          :rich_sequence_training_path => options[:training_sequence_path],
                                          :glimmer_dir_path => options[:glimmer_dir],
                                          :suppress_messages => true)
          end
        else
          if options[:verbose]
            puts "Predicting genes for file #{sequence_file} using iterated glimmer...."
          else
            print "."
          end
          predict_using_iterated_glimmer(:suppress_messages => true, :input_sequence_path => sequence_file, :glimmer_predict_filename => File.basename(sequence_file, File.extname(sequence_file)),:glimmer_dir_path => options[:glimmer_dir])
          predict_file = File.basename(sequence_file, File.extname(sequence_file)) + ".predict"
        end
        if options[:verbose]
          puts "Converting #{sequence_file} glimmer prediction to a genbank file ...."
        else
          print "."
        end
        glimmer_genbank_file = glimmer_prediction_to_rich_sequence_file(:suppress_messages => true, :glimmer_predict_file => predict_file, :input_sequence_path  => sequence_file)
        files_with_cds << glimmer_genbank_file
      else
        files_with_cds << sequence_file
      end
    end

    cds_multi_fasta_protein_file = File.open(options[:cds_multi_fasta_file], "w")
    read_cds_and_write_to_file(files_with_cds, cds_multi_fasta_protein_file)
    processing_indicator(5)

    cds_multi_fasta_protein_file.close
  end
  # method to perform cdhit using heirachical cutoffs
  # @param [Hash] options The options for the method
  # @option options [String] :root_folder the path and name of the folder where the results files will be written
  # @option options [String] :cds_multi_fasta_file The name of the multi fasta that will be created
  # @option options [Array] :clustering_cutoffs The percent id cutoffs to use when performing heirachical clustering with cd-hit
  # @option options [String] :cdhit_dir The directory containing the cd-hit executable
  def make_heirachical_clusters(options) # perform heirachical clustering with cd hit
    default_options = {
      :clustering_cutoffs => [99, 98, 95, 90, 85],
      :cdhit_dir => "/usr/local/cdhit",
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :root_folder, :required => true, :type => :string
      option :cds_multi_fasta_file, :required => true, :type => :string
    end


    Dir.chdir(options[:root_folder])
    options[:clustering_cutoffs].each_with_index do |clustering_cutoff, index|
      puts "Clustering proteins with cutoff #{clustering_cutoff}"
      if index == 0
        `#{options[:cdhit_dir]}/cd-hit -i #{options[:cds_multi_fasta_file]} -o cdhit_clusters-#{clustering_cutoff} -c #{clustering_cutoff/100.to_f} -n 5 -d 200`
      else
        `#{options[:cdhit_dir]}/cd-hit -i cdhit_clusters-#{options[:clustering_cutoffs][index-1]} -o cdhit_clusters-#{clustering_cutoff} -c #{clustering_cutoff/100.to_f} -n 5 -d 200`
      end
    end
  end

  # method to construct a database from cd-hit results
  # @param [Hash] options The options for the method
  # @option options [String] :root_folder the path and name of the folder where the results files will be written
  # @option options [String] :db_location The path and name of the database to be created
  # @option options [Array] :strain_names An array of strain names in the same order as the sequence files supplied to be clustered with cd-hit
  # @option options [Array] :clustering_cutoffs The percent id cutoffs to use when performing heirachical clustering with cd-hit
  # @option options [Boolean] :verbose If true then for debug purposes the time taken to create 100 clusters will be printed to screen
  def make_db_clusters(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'
    require 'core-access/cluster_classes'
    default_options = {
      :clustering_cutoffs => [99, 98, 95, 90, 85],
      :verbose => false
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :strain_names, :required => true, :type => :array
      option :db_location, :required => true, :type => :string
      option :root_folder, :required => true, :type => :string
    end

    Dir.chdir(options[:root_folder])

    connection = make_db_connection(options[:db_location]).connection
    if options[:verbose]
      make_cluster_db_schema # make the datbase tables
    else
      silence do
        make_cluster_db_schema # make the datbase tables
      end
    end
              
    cluster_collection = ClusterCollection.new
    options[:clustering_cutoffs].reverse.each do |cutoff|
      puts "Processing cdhit clustering results with cutoff #{cutoff}"
      # process cdhit file
      clusters_details = read_clstr_file("cdhit_clusters-#{cutoff}.clstr")
      # add clusters to ClusterCollection
      process_clusters(clusters_details, cluster_collection)
    end

    # create strains in database
    options[:strain_names].each do |strain_name|
      Strain.create(:name => strain_name)
    end

    # create clusters
    # create_clusters_in_db_by_sql(connection, cluster_collection, options[:clustering_cutoffs].last, options[:verbose])
    create_clusters_in_db(cluster_collection, options[:clustering_cutoffs].last, options[:verbose])


    puts "updating cluster number_of_members and number_of_strains"
    # update cluster info
    update_cluster_number_of_members(Cluster.find(:all))
    update_cluster_number_of_strains(Cluster.find(:all))

  end
  # method to construct a database from cd-hit results
  # @param [Hash] options The options for the method
  # @option options [String] :root_folder the path and name of the folder where the results files will be written
  # @option options [String] :db_location The path and name of the database to be created
  # @option options [Array] :strain_names An array of strain names in the same order as the sequence files supplied to be clustered with cd-hit
  # @option options [Array] :sequence_files An array of filepaths listing each of the input sequence files. Should be the same order as the strain names
  # @option options [Bolean] :verbose Whether to ouput verbose debug information
  def add_representative_sequences_to_cluster(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'
    require 'bioutils/rich_sequence_utils'

    default_options = {
      :verbose => false
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :db_location, :required => true, :type => :string
      option :root_folder, :required => true, :type => :string
      option :sequence_files, :required => true, :type => :array
      option :strain_names, :required => true, :type => :array
    end

    genbank_files = Array.new
    options[:sequence_files].each do |sequence_file|
      if guess_sequence_format(sequence_file) == :fasta
        genbank_files << "#{options[:root_folder]}/#{File.basename(sequence_file, File.extname(sequence_file)) + "_glimmer_genes.gbk"}"
      else
        genbank_files << sequence_file
      end
    end

    coding_sequences = load_genes_from_genomes(genbank_files, options[:strain_names], options[:verbose])
    connection = make_db_connection(options[:db_location]).connection
    add_representatives_to_cluster(coding_sequences)
  end

  # method to construct clusters of clusters at a specified cutoff (here called super clusters)
  # @param [Hash] options The options for the method
  # @option options [String] :db_location The path and name of the database to be created
  # @option options [Integer] :cutoff The precent cutoff used to cluster the clusters with cd-hit
  # @option options [String] :cdhit_dir The directory containing the cd-hit executable
  def make_super_clusters(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'
    require 'tempfile'

    default_options = {
      :cutoff  => 70,
      :cdhit_dir => "/usr/local/cdhit"
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :db_location, :required => true, :type => :string
    end

    connection = make_db_connection(options[:db_location]).connection

    tmpfile_object = write_cluster_representatives_to_file

    # make cdhit clusters
    `#{options[:cdhit_dir]}/cd-hit -i #{tmpfile_object.path} -o /tmp/cdhit_super_clusters-#{options[:cutoff]} -c #{options[:cutoff]/100.to_f} -n 5  -d 200`

    # read in cluster info from cdhit output
    super_clusters = read_super_clstr_file("/tmp/cdhit_super_clusters-#{options[:cutoff]}.clstr")
    add_super_clusters_to_db(super_clusters, options[:cutoff])

    puts "updating cluster number_of_members and number_of_strains"
    # update cluster info
    update_cluster_number_of_members(Cluster.where("clusters.is_parent_cluster = ?", true).all)
    update_cluster_number_of_strains(Cluster.where("clusters.is_parent_cluster = ?", true).all)

    tmpfile_object.close
  end

  ################################ sub methods #############################

  def read_cds_and_write_to_file(sequence_files, multi_fasta_ouput_file)
    require 'bioutils/rich_sequence_utils'
    genome_counter = 0
    sequence_files.each do |sequence_file|
      genome_counter += 1
      cds_counter = 0

      sequence_format = guess_sequence_format(sequence_file)

      case sequence_format
      when :genbank
        sequence_flatfile = Bio::FlatFile.open(Bio::GenBank,sequence_file)
      when :embl
        sequence_flatfile = Bio::FlatFile.open(Bio::EMBL,sequence_file)
      else
        puts "All sequence files should be of genbank or embl format"
        exit
      end
      while sequence_entry = sequence_flatfile.next_entry
        sequence = sequence_entry.seq
        sequence_entry.each_cds do |cds|
          cds_counter += 1
          if cds_counter % 500 == 0
            processing_indicator(cds_counter/500 % 4)
          end
          location_string =  cds.locations.to_s
          dna_cds_sequence  = sequence.splice(location_string)
          protein_cds_sequence = dna_cds_sequence.translate
          unless dna_cds_sequence =~ /nnnnnnnnnn/i # ignore genes with N's (e.g from scaffolds)
            multi_fasta_ouput_file.puts ">#{genome_counter}-#{cds_counter}"
            multi_fasta_ouput_file.puts protein_cds_sequence
          end
        end
      end
    end
  end

  def read_clstr_file(clstr_file_path)#read clstr output file and record cluster members
    clstr_file_handle = File.open(clstr_file_path)
    clusters = Array.new
    cluster = nil

    while line = clstr_file_handle.gets
      line.chomp!
      if line =~ />Cluster/ # new cluster
        unless cluster.nil?
          clusters << cluster
        end
        cluster = Array.new
      else
        match_data = line.match(/\d+?\s+?(\d+)aa, >(\d+?-\d+?)\.\.\.\s(.+)/)
        size = match_data.captures[0].to_i
        name = match_data.captures[1]
        clustering_status = match_data.captures[2]
        if clustering_status == "*"
          clustering_status = "founder"
        else
          clustering_status =~ /at\s((\d|\.)+)%/
          clustering_status = $1.to_f
        end
        cluster << {:name => name, :size => size, :clustering_status => clustering_status}
      end
    end
    clusters << cluster # add final cluster
    return clusters # clusters is a format [[{:name => name1, :size => size1, :clustering_status => clustering_status1}, {:name => name2, :size => size2, :clustering_status => clustering_status2}], [{:name => name3, :size => size3, :clustering_status => clustering_status3}]] An array of arrays where each array member is an array of hash object representing cluster members with :name and :size keys
  end

  def process_clusters(clusters_details, cluster_collection) # process clusters
    if cluster_collection.clusters.keys.empty?
      next_index = 1
    else
      next_index = cluster_collection.clusters.keys.sort.last + 1
    end
    
    clusters_details.each do |cluster_detail| # loop through clusters
      cluster_membership = Array.new
      cluster_detail.each do |cluster_member_detail| # loop through cluster members and check which are already in a cluster
        existing_member_details = cluster_collection.cluster_members[cluster_member_detail[:name]]
        if existing_member_details.nil?
          cluster_index = nil
          member = nil
        else
          cluster_index = existing_member_details[:cluster_index]
          member = existing_member_details[:member]
        end
        
        if member.nil?
          cluster_membership << {:cluster_index => cluster_index, :member => ClusterMember.new(cluster_member_detail[:name], cluster_member_detail[:size], cluster_member_detail[:clustering_status])}
        else
          cluster_membership << {:cluster_index => cluster_index, :member => member}
        end
      end
      
      if cluster_membership.map{|cm| cm[:cluster_index]}.compact.empty? # no existing clusters so create a new one
        new_cluster = CDHitCluster.new
        cluster_collection.clusters[next_index] = new_cluster
        
        cluster_membership.each do |cm|
          new_cluster.add_member(cm[:member])
          cluster_collection.cluster_members[cm[:member].name] = {:member => cm[:member], :cluster_index => next_index}
        end
        
        next_index += 1
      else
        existing_cluster_indicices = cluster_membership.map{|cm| cm[:cluster_index]}.compact.uniq
        merged_cluster_index = existing_cluster_indicices.first
        merged_cluster = cluster_collection.clusters[merged_cluster_index]
        cluster_membership.each do |cm|
          if cm[:cluster_index] == merged_cluster_index
            # do nothing member is already in merged cluster
          else# new cluster member or from existing cluster that will be deleted once merged
            merged_cluster.add_member(cm[:member])
            cluster_collection.cluster_members[cm[:member].name] = {:member => cm[:member], :cluster_index => merged_cluster_index}
          end
        end
        existing_cluster_indicices.slice(1..-1).each do |cluster_index_to_delete| # delete other clusters
          cluster_collection.clusters.delete(cluster_index_to_delete)
        end
      end
    end
  end

  def create_clusters_in_db(cluster_collection, cutoff, show_timings = false)
    if show_timings
      cluster_counter = 0
      a = Time.now
    end
    cluster_collection.clusters.each do |cluster_index, cluster|
      ActiveRecord::Base.transaction do
        if show_timings
          cluster_counter += 1
          if cluster_counter % 100 == 0
            b = Time.now
            mm, ss = (b-a).divmod(60)
            puts "Adding 100 clusters took #{mm} minutes #{ss} seconds"
            a = Time.now
          end
        end
        # puts cluster_index
        db_cluster = Cluster.new
        db_cluster.cutoff = cutoff
        db_cluster.is_parent_cluster = false
        db_cluster.save
        # genes = Array.new
        previous_strain_number = 0
        previous_gene_number = 0
        gene = nil
        cluster.members.sort_by{|m| [m.name.split("-").first.to_i, m.name.split("-").last.to_i]}.each do |member|
          strain_number, gene_number = member.name.split("-").map{|e| e.to_i}
          if strain_number == previous_strain_number && gene_number == previous_gene_number + 1 # should probably be joined
           gene.name = gene.name.to_s + "-#{gene_number}"
           gene.aa_length += member.size
           gene.save
          else
            gene = Gene.new
            gene.name = gene_number
            gene.strain = Strain.find(strain_number)
            gene.aa_length = member.size
            gene.save
            ClusterMembership.create(:gene_id => gene.id, :cluster_id => db_cluster.id, :status => member.clustering_status)
            # genes << gene
          end
          previous_strain_number = strain_number
          previous_gene_number = gene_number
        end
        # db_cluster.genes = genes
        # db_cluster.save
        puts "Cluster #{db_cluster.id} created" if db_cluster.id % 100 == 0
      end
    end
  end

  def create_clusters_in_db_by_sql(connection, cluster_collection, cutoff, show_timings = false)
    if show_timings
      cluster_counter = 0
      a = Time.now
    end
    cluster_collection.clusters.each do |cluster_index, cluster|
      if show_timings
        cluster_counter += 1
        if cluster_counter % 100 == 0
          b = Time.now
          mm, ss = (b-a).divmod(60)
          puts "Adding 100 clusters took #{mm} minutes #{ss} seconds"
          a = Time.now
        end
      end
      # puts cluster_index
      connection.execute("INSERT INTO clusters (cutoff, is_parent_cluster) VALUES (#{cutoff}, 'f')")
      result = connection.execute("SELECT * FROM clusters ORDER BY id DESC LIMIT 1").first
      cluster_id = result["id"]

      previous_strain_number = 0
      previous_gene_number = 0
      gene_id = 0
      cluster.members.sort_by{|m| [m.name.split("-").first.to_i, m.name.split("-").last.to_i]}.each do |member|
        strain_number, gene_number = member.name.split("-").map{|e| e.to_i}
        if strain_number == previous_strain_number && gene_number == previous_gene_number + 1 # should probably be joined
          result = connection.execute("SELECT * FROM genes WHERE id = #{gene_id}").first
          gene_name = result["name"]
          connection.execute("UPDATE genes SET name = '#{gene_name + "-" + gene_number}' WHERE id = #{gene_id}")
        else
          connection.execute("INSERT INTO genes (name, strain_id) VALUES ('#{gene_number}', #{strain_number})")
          result = connection.execute("SELECT * FROM genes ORDER BY id DESC LIMIT 1").first
          gene_id = result["id"]
          # join gene to cluster
          connection.execute("INSERT INTO cluster_memberships (cluster_id, gene_id) VALUES (#{cluster_id}, #{gene_id})")
        end
        previous_strain_number = strain_number
        previous_gene_number = gene_number
      end
      puts "Cluster #{cluster_id} created by SQL"
    end
    puts "Committing transaction to database"
  end

  def load_genes_from_genomes(genbank_files, strain_names, verbose = false) # pass in a list of genbank file paths
    coding_sequences = Hash.new
    # load all genes into a hash
    genome_counter = 0
    genbank_files.each do |genbank_file|
      genome_counter += 1
      cds_counter = 0
      genbank_flatfile = Bio::FlatFile.open(Bio::GenBank,genbank_file)
      genome_name = strain_names[genome_counter - 1]
      coding_sequences[genome_name] = Array.new
      offset = 0
      while genbank_entry = genbank_flatfile.next_entry
        if verbose
          puts "Adding genes from #{genbank_entry.definition}"
        end
        genbank_sequence = genbank_entry.seq
        genbank_entry.each_cds do |cds|
          location_string =  cds.locations.to_s
          offset_location_string = location_string.gsub(/(\d+)/){|match| match.to_i + offset}
          dna_cds_sequence  = genbank_sequence.splice(location_string)
          coding_sequences[genome_name] << {:location_string => location_string, :offset_location_string => offset_location_string, :dna_cds_sequence =>  dna_cds_sequence}
        end
        offset += genbank_sequence.length
      end
    end
    return coding_sequences
  end

  def add_representatives_to_cluster(coding_sequences)
    # add representative gene to Clusters
    Cluster.find(:all).each do |cluster|
      ActiveRecord::Base.transaction do
        puts "adding repesentative to #{cluster.id}" if cluster.id % 100 == 0
        longest_gene = nil
        longest_gene_length = 0
        cluster.genes.each do |gene|
          gene_names = gene.name.split("-")
          strain_name = gene.strain.name
          locations = Array.new
          offset_locations = Array.new
          gene_names.each do |gene_name|
            offset_locations << coding_sequences[strain_name][gene_name.to_i - 1][:offset_location_string]
            locations << coding_sequences[strain_name][gene_name.to_i - 1][:location_string]
          end
          gene.location = offset_locations.join(", ")
          gene.relative_location = locations.join(", ") # relative location is the location within the original contig
          gene.save
          length = 0
          locations.each do |location|
            match =location.match(/(\d+)\D*?\.\.\D*?(\d+)/)
            start_pos = match.captures[0].to_i
            end_pos = match.captures[1].to_i
            length += end_pos - start_pos
          end
          if length > longest_gene_length
            longest_gene = gene
          end
        end
        gene_names = longest_gene.name.split("-")
        strain_name = longest_gene.strain.name
        longest_gene.sequence = ""
        gene_names.each do |gene_name|
          longest_gene.sequence += coding_sequences[strain_name][gene_name.to_i - 1][:dna_cds_sequence] # += is to cater for genes that are joined
        end
        longest_gene.save
        cluster.representative = longest_gene
        cluster.save
      end
    end
  end

  def write_cluster_representatives_to_file
    tmpfile_object = Tempfile.new('temp')
    Cluster.all.each do |cluster|
      tmpfile_object.puts ">#{cluster.id}"
      biosequence = Bio::Sequence.new(cluster.representative.sequence)
      biosequence.na
      tmpfile_object.puts biosequence.translate(1,11)
    end
    tmpfile_object.close
    return tmpfile_object
  end

  def read_super_clstr_file(clstr_filepath)
    super_clusters = Array.new
    super_cluster = nil

    clstr_file = File.open(clstr_filepath)

    while line = clstr_file.gets
      line.chomp!
      if line =~ />Cluster/ # new cluster
        unless super_cluster.nil?
          super_clusters << super_cluster
        end
        super_cluster = Array.new
      else
        match_data = line.match(/\d+?\s+?(\d+)aa, >(\d+?)\.\.\.\s(.+)/)
        size = match_data.captures[0].to_i
        name = match_data.captures[1]
        clustering_status = match_data.captures[2]
        if clustering_status == "*"
          clustering_status = "founder"
        else
          clustering_status =~ /at\s((\d|\.)+)%/
          clustering_status = $1.to_f
        end
        super_cluster << {:name => name, :size => size, :clustering_status => clustering_status}
      end
    end
    clstr_file.close
    return super_clusters
  end

  def add_super_clusters_to_db(super_clusters, cutoff)
    # process clusters to make parent child clusters in database
    parent_cluster_name = 0
    super_clusters.each do |super_cluster|
      if super_cluster.size > 1
        parent_cluster_name += 1
        parent_cluster = Cluster.new
        parent_cluster.name = parent_cluster_name
        parent_cluster.cutoff = cutoff
        parent_cluster.is_parent_cluster = true
        parent_cluster.save

        child_cluster_name = 0
        longest_gene_size = 0
        parent_cluster_represenative = nil
        super_cluster.each do |cluster_member|
          child_cluster_name += 1
          child_cluster = Cluster.find(cluster_member[:name].to_i)
          child_cluster.parent_cluster = parent_cluster
          child_cluster.name = "#{parent_cluster_name}.#{child_cluster_name}"
          child_cluster.save

          if child_cluster.representative.sequence.size > longest_gene_size
            parent_cluster_represenative = child_cluster.representative
            longest_gene_size = child_cluster.representative.sequence.size
          end

          child_cluster.genes.each do |gene|
            ClusterMembership.create(:gene_id => gene.id, :cluster_id => parent_cluster.id, :status => cluster_member[:clustering_status])
          end
        end
        parent_cluster.representative = parent_cluster_represenative
        parent_cluster.save
        puts "Created new parent cluster #{parent_cluster.id}"

        # puts super_cluster.map{|c| c[:name]}.join(" ")
      end
    end
  end

  def add_super_clusters_to_db_by_sql(connection, super_clusters, cutoff)
    # process clusters to make parent child clusters in database
    parent_cluster_name = 0
    super_clusters.each do |super_cluster|
      if super_cluster.size > 1
        parent_cluster_name += 1
        connection.execute("INSERT INTO clusters (name, cutoff, is_parent_cluster) VALUES ('#{parent_cluster_name}', #{cutoff}, 't')")
        result = connection.execute("SELECT * FROM clusters ORDER BY id DESC LIMIT 1").first
        parent_cluster_id = result["id"]

        child_cluster_name = 0
        longest_gene_size = 0
        parent_cluster_represenative_id = nil
        parent_cluster_gene_ids = Array.new
        super_cluster.each do |cluster_member|
          child_cluster_name += 1
          connection.execute("UPDATE clusters SET name = '#{parent_cluster_name}.#{child_cluster_name}', parent_id = #{parent_cluster_id} WHERE id = #{cluster_member[:name]}")
          child_cluster = Cluster.find(cluster_member[:name].to_i)

          result = connection.execute("SELECT * FROM clusters WHERE id = #{cluster_member[:name]}").first
          representative_id = result["representative_id"]
          result = connection.execute("SELECT * FROM genes WHERE id = #{representative_id}").first
          representative_sequence = result["sequence"]
          if representative_sequence.size > longest_gene_size
            parent_cluster_represenative_id = representative_id
            longest_gene_size = representative_sequence.size
          end

          result = connection.execute("SELECT * FROM cluster_memberships WHERE cluster_id = #{cluster_member[:name]}")
          result.each do |cluster_membership|
            parent_cluster_gene_ids << cluster_membership["gene_id"]
          end
        end
        # add representative
        connection.execute("UPDATE clusters SET representative_id = #{parent_cluster_represenative_id} WHERE id = #{parent_cluster_id}")
        # add genes to parent cluster
        parent_cluster_gene_ids.each do |parent_cluster_gene_id|
          connection.execute("INSERT INTO cluster_memberships (cluster_id, gene_id) VALUES (#{parent_cluster_id}, #{parent_cluster_gene_id})")
        end
        puts "Created new parent cluster #{parent_cluster_id}"

        # puts super_cluster.map{|c| c[:name]}.join(" ")
      end
    end
  end

  def update_cluster_number_of_members(clusters)
    clusters.each do |cluster|
      cluster.number_of_members = cluster.genes.size
      cluster.save
    end
  end

  def update_cluster_number_of_strains(clusters)
    clusters.each do |cluster|
      cluster.number_of_strains = Cluster.number_of_strains(cluster.id)
      cluster.save
    end
  end

end
