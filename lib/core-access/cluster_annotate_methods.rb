module ClusterAnnotate
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
      :local_db_blast_program => "blastn",
      :local_db_percent_identity_cutoff => 85,
      :local_db_minimum_hit_length => 85,
      :ncbi_blast_program => "blastp",
      :ncbi_percent_identity_cutoff => 80,
      :ncbi_minimum_hit_length => 80,
      :annotate_vs_local_db => true,
      :annotate_by_remote_blast => true
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
    Cluster.all.each do |cluster|
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
          puts "Annotating cluster #{cluster_id}"
          options[:representative_biosequence] = representative_biosequence
          reciprocal_hit_details = annotate_cluster_sequence(options)
          cluster_annotation_array << [cluster_id, reciprocal_hit_details]
        end
        cluster_annotation_array
      end
      # cleanup collected annotation
      collected_annotations.each do |collected_annotation|\
        collected_annotation.each do |cluster_annotation|
          cluster_annotations[cluster_annotation.first] = cluster_annotation.last # cluster_id as key , annotations as value
        end
      end
    else
      # serial annotation
      representative_sequences.each do |cluster_id, representative_biosequence|
        puts "Annotating cluster #{cluster_id}"
        options[:representative_biosequence] = representative_biosequence
        reciprocal_hit_details = annotate_cluster_sequence(options)
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

  def annotate_cluster_sequence(options)
    reciprocal_hit_details = GenomeReciprocalHitAnnotator.annotate_sequence(
                :biosequence  => options[:representative_biosequence],
                :query_sequence_database_name => "blast_databases/cluster_representatives",
                :accept_reciprocal_hits_with_multiple_hits_incl_query => options[:accept_reciprocal_hits_with_multiple_hits_incl_query],
                :accept_first_reciprocal_hit_containing_query => options[:accept_first_reciprocal_hit_containing_query],
                :blast_dir  => options[:blast_dir],
                :fastacmd_dir => options[:fastacmd_dir],
                :reference_sequences_database_name => "blast_databases/reference_genomes",
                :reference_blast_program => options[:reference_blast_program],
                :reference_percent_identity_cutoff => options[:reference_percent_identity_cutoff],
                :reference_minimum_hit_length => options[:reference_minimum_hit_length],
                :local_blast_db => options[:local_blast_db],
                :local_db_blast_program => options[:local_db_blast_program],
                :local_db_percent_identity_cutoff => options[:local_db_percent_identity_cutoff],
                :local_db_minimum_hit_length => options[:local_db_minimum_hit_length],
                :ncbi_blast_program => options[:ncbi_blast_program],
                :ncbi_percent_identity_cutoff => options[:ncbi_percent_identity_cutoff],
                :ncbi_minimum_hit_length => options[:ncbi_minimum_hit_length],
                :annotate_vs_local_db => options[:annotate_vs_local_db],
                :annotate_by_remote_blast => options[:annotate_by_remote_blast])
    reciprocal_hit_details.delete_if{|hd| hd[0] =~ /(translation|transl_table)/} unless reciprocal_hit_details.nil?
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
end