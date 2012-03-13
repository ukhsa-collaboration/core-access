module Clustering
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

  def create_cds_multi_fasta_file(options)
    require 'bioutils/rich_sequence_utils'
    require 'bioutils/glimmer'
    extend Glimmer

    default_options = {
      :cds_multi_fasta_file => "cds_proteins.fas"
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
          run_glimmer_using_model(:input_sequence_path  => sequence_file, :prefix => options[:training_model_prefix],:suppress_messages => true)
          predict_file = File.basename(sequence_file, File.extname(sequence_file)) + "_glimmer.predict"
        elsif options[:training_sequence_path]
          model_file_prefix = File.basename(options[:training_sequence_path], File.extname(options[:training_sequence_path])) + "_glimmer"
          if File.exists?(model_file_prefix + ".icm")
            puts "Predicting genes for file #{sequence_file} using training model ...."
            run_glimmer_using_model(:input_sequence_path  => sequence_file, :prefix => model_file_prefix,:suppress_messages => true)
            predict_file = File.basename(sequence_file, File.extname(sequence_file)) + "_glimmer.predict"
          else
            puts "Predicting genes for file #{sequence_file} using training sequence ...."
            predict_file = predict_genes_using_glimmer(:input_sequence_path  => sequence_file,
                                          :rich_sequence_training_path => options[:training_sequence_path],
                                          :suppress_messages => true)
          end
        else
          puts "Predicting genes for file #{sequence_file} using iterated glimmer...."
          predict_using_iterated_glimmer(:suppress_messages => true, :input_sequence_path => sequence_file, :glimmer_predict_filename => File.basename(sequence_file, File.extname(sequence_file)))
          predict_file = File.basename(sequence_file, File.extname(sequence_file)) + ".predict"
        end
        puts "Converting #{sequence_file} glimmer prediction to a genbank file ...."
        glimmer_genbank_file = glimmer_prediction_to_rich_sequence_file(:suppress_messages => true, :glimmer_predict_file => predict_file, :input_sequence_path  => sequence_file)
        files_with_cds << glimmer_genbank_file
      else
        files_with_cds << sequence_file
      end
    end

    cds_multi_fasta_protein_file = File.open(options[:cds_multi_fasta_file], "w")
    read_cds_and_write_to_file(files_with_cds, cds_multi_fasta_protein_file)

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
  # @option options [Boolean] :show_timings If true then for debug purposes the time taken to create 100 clusters will be printed to screen
  def make_db_clusters(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'
    require 'core-access/cluster_classes'
    default_options = {
      :clustering_cutoffs => [99, 98, 95, 90, 85],
      :show_timings => false
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :strain_names, :required => true, :type => :array
      option :db_location, :required => true, :type => :string
      option :root_folder, :required => true, :type => :string
    end

    Dir.chdir(options[:root_folder])

    connection = make_db_connection(options[:db_location]).connection
    make_cluster_db_schema # make the datbase tables
              
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
    # create_clusters_in_db_by_sql(connection, cluster_collection, options[:clustering_cutoffs].last, options[:show_timings])
    create_clusters_in_db(cluster_collection, options[:clustering_cutoffs].last, options[:show_timings])


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
  def add_representative_sequences_to_cluster(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'
    require 'bioutils/rich_sequence_utils'

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

    coding_sequences = load_genes_from_genomes(genbank_files, options[:strain_names])
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

    Dir.chdir("/tmp")

    connection = make_db_connection(options[:db_location]).connection

    tmpfile_object = write_cluster_representatives_to_file

    # make cdhit clusters
    `#{options[:cdhit_dir]}/cd-hit -i #{tmpfile_object.path} -o cdhit_super_clusters-#{options[:cutoff]} -c #{options[:cutoff]/100.to_f} -n 5  -d 200`

    # read in cluster info from cdhit output
    super_clusters = read_super_clstr_file("cdhit_super_clusters-#{options[:cutoff]}.clstr")
    add_super_clusters_to_db(super_clusters, options[:cutoff])

    puts "updating cluster number_of_members and number_of_strains"
    # update cluster info
    update_cluster_number_of_members(Cluster.where("clusters.is_parent_cluster = ?", true).all)
    update_cluster_number_of_strains(Cluster.where("clusters.is_parent_cluster = ?", true).all)

    tmpfile_object.close
  end

  def annotate_clusters(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'
    require 'genome/genome_reciprocal_hit_annotator'
    require 'bioutils/blast'
    extend Blast

    default_options = {
      :reference_blast_program => "blastn",
      :reference_percent_identity_cutoff => 95,
      :reference_minimum_hit_length => 85,
      :microbial_genomes_blast_program => "blastn",
      :microbial_genomes_percent_identity_cutoff => 85,
      :microbial_genomes_minimum_hit_length => 85,
      :ncbi_blast_program => "blastp",
      :ncbi_percent_identity_cutoff => 80,
      :ncbi_minimum_hit_length => 80
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :db_location, :required => true, :type => :string
      option :root_folder, :required => true, :type => :string
    end

    Dir.chdir(options[:root_folder])

    connection = make_db_connection(options[:db_location]).connection
    unless options[:reference_database_path]
      puts "Making reference sequence blast databases"
      # make protein blast reference database
       blast_database_from_rich_sequences(:input_sequences => options[:reference_genomes], :database_labels => :full_annotation, :database_name => "reference_genomes", :protein => true, :final_db_location => "blast_databases", :formatdb_dir => options[:blast_dir])
      # make nucleotide blast reference database
      blast_database_from_rich_sequences(:input_sequences => options[:reference_genomes], :database_labels => :full_annotation, :database_name => "reference_genomes", :protein => false, :final_db_location => "blast_databases", :formatdb_dir => options[:blast_dir])
      options[:reference_database_path] = "blast_databases/reference_genomes"
    end
    puts "Making blast database from cluster reference sequences"
    make_blast_databases_from_clusters

    representative_sequences = Array.new
    Cluster.order("id DESC").all.each do |cluster|
      next if cluster.representative.nil? # skip those clusters without a repesentative (super clusters)
      next unless cluster.representative.annotations.empty? # skip those clusters already annotated
      representative_cds = cluster.representative
      representative_biosequence = Bio::Sequence.new(representative_cds.sequence)
      representative_biosequence.entry_id = cluster.id.to_s
      representative_biosequence.na
      representative_sequences << [cluster.id, representative_biosequence]
    end

    cluster_annotations = Hash.new

    # perform annotation
    if options[:parallel_processors]
      # parallel annotation
      require 'forkoff'
      representative_sequences_slices = Array.new
      representative_sequences.each_slice(10) do |representative_sequences_slice|
        representative_sequences_slices << representative_sequences_slice
      end
      # parallel loop starts here
      collected_annotations = representative_sequences_slices.forkoff :processes => options[:parallel_processors], :strategy => :file do |*representative_sequences_slice|
        cluster_annotation_array = Array.new
        representative_sequences_slice.each do |cluster_id, representative_biosequence|
          reciprocal_hit_details = annotate_cluster_sequence(cluster_id, representative_biosequence)
          cluster_annotation_array << [cluster_id, reciprocal_hit_details]
        end
        cluster_annotation_array
      end
      # cleanup collected annotation
      collected_annotations.each do |collected_annotation|
        collected_annotation.each do |cluster_annotation|
          cluster_annotations[cluster_annotation.first] = cluster_annotation.last # cluster_id as key , annotations as value
        end
      end
    else
      # serial annotation
      representative_sequences.each do |cluster_id, representative_biosequence|
        reciprocal_hit_details = annotate_cluster_sequence(cluster_id, representative_biosequence)
        cluster_annotations[cluster_id] = reciprocal_hit_details # cluster_id as key , annotations as value
      end

    end
   


    # apply annotations to database
    cluster_annotations.each do |cluster_id, cluster_annotation|
      next if cluster_annotation.nil?
      cluster = Cluster.find(cluster_id)
      representative_cds = cluster.representative
      cluster_annotation.each do |annotation|
        representative_cds.annotations << Annotation.new(:qualifier => annotation.first, :value => annotation.last )
      end
    end
  end

  def find_shared_clusters(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'

    default_options = {
      :unique  => false,
      :excluded_cluster_ids => false,
      :include_parent_clusters => false
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :strain_names, :required => true, :type => :array
    end

    test_connection(options[:db_location])

    *strain_names = options[:strain_names]
    where_statement = strain_names.collect{|strain_name| "strains.name = '#{strain_name}' OR "}.join("").sub(/ OR $/, "")

    where_statement += " AND clusters.id NOT IN (#{options[:excluded_cluster_ids].join(", ")})" if options[:excluded_cluster_ids]

    where_statement += " AND clusters.is_parent_cluster = 'f'" unless options[:include_parent_clusters]

    sql_statement = "SELECT * FROM (SELECT clusters.* FROM clusters INNER JOIN cluster_memberships ON clusters.id = cluster_memberships.cluster_id INNER JOIN genes ON genes.id = cluster_memberships.gene_id INNER JOIN strains ON strains.id = genes.strain_id WHERE (#{where_statement}) GROUP BY clusters.id, strain_id)"

    sql_statement += " WHERE number_of_members = #{strain_names.size}" if options[:unique]
    sql_statement += " GROUP BY id HAVING COUNT(*) = #{strain_names.size}"

    return Cluster.find_by_sql(sql_statement)
  end

  def find_unique_clusters(options)
    options.merge!(:unique => true)
  end

  def find_core_clusters(options = {})
    default_options = {
      :db_location => nil
    }
    options.reverse_merge!(default_options)
    
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'

    test_connection(options[:db_location])

    strain_names = Strain.all.map{|strain| strain.name}
    find_shared_clusters(options.merge(:strain_names => strain_names))
  end

  def find_genes_with_annotation(options)
    options = MethodArgumentParser::Parser.check_options options  do
      option :annotation, :required => true, :type => :string
    end

    test_connection

    Gene.joins(:annotations).where("annotations.value LIKE ?", "%#{options[:annotation]}%")
  end

  def find_clusters_having_genes(options)
    options = MethodArgumentParser::Parser.check_options options  do
      option :gene_ids, :required => true, :type => :array
    end

    Cluster.joins(:genes).where("genes.id IN (?)", options[:gene_ids]).select("DISTINCT clusters.*")
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

    # ActiveRecord::Base.transaction do
      # Remove transactions
      # ActiveRecord::ConnectionAdapters::SQLiteAdapter.class_eval do
      #   def begin_db_transaction
      #   end

      #   def commit_db_transaction
      #   end
      # end
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
      puts "Committing transaction to database"
    # end
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

  def load_genes_from_genomes(genbank_files, strain_names) # pass in a list of genbank file paths
    coding_sequences = Hash.new
    # load all genes into a hash
    genome_counter = 0
    genbank_files.each do |genbank_file|
      genome_counter += 1
      cds_counter = 0
      genbank_flatfile = Bio::FlatFile.open(Bio::GenBank,genbank_file)
      genome_name = strain_names[genome_counter - 1]
      coding_sequences[genome_name] = Array.new
      while genbank_entry = genbank_flatfile.next_entry
        puts "Adding genes from #{genbank_entry.definition}"
        genbank_sequence = genbank_entry.seq
        genbank_entry.each_cds do |cds|
          location_string =  cds.locations.to_s
          dna_cds_sequence  = genbank_sequence.splice(location_string)
          coding_sequences[genome_name] << [location_string, dna_cds_sequence]
        end
      end
    end
    return coding_sequences
  end

  def add_representatives_to_cluster(coding_sequences)
    # add representative gene to Clusters
    Cluster.find(:all).each do |cluster|
      puts "adding repesentative to #{cluster.id}" if cluster.id % 100 == 0
      longest_gene = nil
      longest_gene_length = 0
      cluster.genes.each do |gene|
        gene_names = gene.name.split("-")
        strain_name = gene.strain.name
        locations = Array.new
        gene_names.each do |gene_name|
          locations << coding_sequences[strain_name][gene_name.to_i - 1].first
        end
        # puts locations.join(", ")
        gene.location = locations.join(", ")
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
        longest_gene.sequence += coding_sequences[strain_name][gene_name.to_i - 1].last # += is to cater for genes that are joined
      end
      longest_gene.save
      cluster.representative = longest_gene
      cluster.save
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

  def make_blast_databases_from_clusters
    require 'tempfile'
    require 'bioutils/blast'
    extend Blast
    nucleotide_query_sequences = Tempfile.new('temp')
    protein_query_sequences = Tempfile.new('temp')
    
    Cluster.all.each do |cluster|
      next if cluster.representative.nil? # skip those clusters without a repesentative (super clusters)
      representative_cds = cluster.representative
      biosequence = Bio::Sequence.new(representative_cds.sequence)
      biosequence.na
      nucleotide_query_sequences.puts ">#{cluster.id}"
      nucleotide_query_sequences.puts biosequence
      
      protein_query_sequences.puts ">#{cluster.id}"
      protein_query_sequences.puts biosequence.translate(1,11)
    end
    
    nucleotide_query_sequences.close
    protein_query_sequences.close
    
    formatdb(:path_to_fasta_input_sequence => nucleotide_query_sequences.path, :database_name => "cluster_representatives", :final_db_location => "blast_databases")
    
    formatdb(:path_to_fasta_input_sequence => protein_query_sequences.path, :database_name => "cluster_representatives", :formatdb_options => "-o T -p T", :final_db_location => "blast_databases")
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

  def annotate_cluster_sequence(cluster)
    puts "Annotating cluster #{cluster_id}"
    reciprocal_hit_details = GenomeReciprocalHitAnnotator.annotate_sequence(
                :biosequence  => representative_biosequence,
                :reference_blast_program => options[:reference_blast_program],
                :reference_percent_identity_cutoff => options[:reference_percent_identity_cutoff],
                :reference_minimum_hit_length => options[:reference_minimum_hit_length],
                :local_microbial_blast_DB => options[:microbial_genomes_blast_db],
                :microbial_genomes_blast_program => options[:microbial_genomes_blast_program],
                :microbial_genomes_percent_identity_cutoff => options[:microbial_genomes_percent_identity_cutoff],
                :microbial_genomes_minimum_hit_length => options[:microbial_genomes_minimum_hit_length],
                :ncbi_blast_program => options[:ncbi_blast_program],
                :ncbi_percent_identity_cutoff => options[:ncbi_percent_identity_cutoff],
                :ncbi_minimum_hit_length => options[:ncbi_minimum_hit_length])
    reciprocal_hit_details.delete_if{|hd| hd[0] =~ /(translation|transl_table)/} unless reciprocal_hit_details.nil?

end
