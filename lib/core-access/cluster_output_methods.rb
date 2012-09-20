module ClusterOutput
  # a method to print summary information about a cluster
  # @param [Cluster] cluster The cluster object (ActiveRecord) to return summary information for
  # @param [Array] attributes An array of strings listing the cluster attributes that should be returned in the summary
  # default is ["id", "cutoff", "is_parent_cluster", "number_of_members", "number_of_strains"])
  def print_cluster_summary(cluster, attributes = ["id", "cutoff", "is_parent_cluster", "number_of_members", "number_of_strains"])
    cluster_summary_array = Array.new
    attributes.each do |field|
      cluster_summary_array << cluster.send(field)
    end
    puts cluster_summary_array.join("\t")
  end

  # a method to output presence and absence data from clusters as a tab-separated 1/0 matrix with each
  # row representing a cluster and each column the presence of a gene from strain x in that cluster.
  # The first item in each row will be the cluster id and product if the cluster has such an annotation
  # @param [Hash] options A hash of options
  # @option options [String] :db_location The path and name of the core-access database
  # @option options [String] :output_filepath The path of the output file that wiil be created
  # containing presence and absence data 
  # @option options [Boolean] :without_core_genes Whether to include core genes in the output. Generally this
  # will be inadvisable since every core gene will consist of a row of all 1's. To find core genes use the 
  # find_core_clusters instead
  def output_gene_presence_absence(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'

    default_options = {
      :without_core_genes  => true,
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :db_location, :required => true, :type => :string
      option :output_filepath, :required => true, :type => :string
    end
    
    test_connection(options[:db_location])

    output_file = File.open(options[:output_filepath], "w")

    strains = Strain.all
    if options[:annotation_order]
      output_file.print "Cluster id\t"
      options[:annotation_order].each do |annotation_qualifier|
        output_file.print "#{annotation_qualifier}\t"
      end
      output_file.puts strains.map{|strain| strain.name}.join("\t")
    else
      output_file.puts "Descriptor (id: gene; protein_id; product)\t#{strains.map{|strain| strain.name}.join("\t")}"
    end

    clusters = get_clusters(options)

    counter = 0
    clusters.each do |cluster|
      ActiveRecord::Base.transaction do
        counter += 1
        puts "Completed output for #{counter} clusters" if counter % 100 == 0
        if options[:annotation_order]
          cluster_descriptor = "#{cluster.id.to_s}"
          options[:annotation_order].each do |annotation_qualifier|
            cluster_descriptor += "\t"
            annotations = cluster.representative.annotations.select{|annotation| annotation.qualifier == annotation_qualifier}
            cluster_descriptor += annotations.map{|annotation| annotation.value}.join(", ") unless annotations.empty?
          end
        else  
          products = cluster.representative.annotations.select{|annotation| annotation.qualifier == "product"}
          genes = cluster.representative.annotations.select{|annotation| annotation.qualifier == "gene"}
          protein_ids = cluster.representative.annotations.select{|annotation| annotation.qualifier == "protein_id"}
          cluster_descriptor = "#{cluster.id.to_s}: "
          if products.empty? && genes.empty? && protein_ids.empty?
            cluster_descriptor += "No annotation"
          else
            cluster_descriptor += genes.map{|gene| gene.value}.join(", ") + "; " unless genes.empty?
            cluster_descriptor += protein_ids.map{|protein_id| protein_id.value}.join(", ") + "; " unless protein_ids.empty?
            cluster_descriptor += products.map{|product| product.value}.join(", ") + "; " unless products.empty?
          end
        end
        output_file.print "#{cluster_descriptor}"
        strains_with_member_in_cluster = Strain.joins(:genes => :clusters).where("clusters.id = ?", cluster.id).group(:id).all
        strains.each do |strain|
          if strains_with_member_in_cluster.include?(strain)
            output_file.print "\t1" #gene presence
          else
            output_file.print "\t0" #gene absence
          end
        end
        output_file.puts
      end
    end
    output_file.close
  end
  
  # This method outputs the gene names for each strain given a list of clusters or if not if no list is given
  # the clusters ordered by id. Since the gene names are sequential numbers based on the order in the original
  # sequence file used to generate the database the output can be used to determine if the order of the genes
  # in a particular sequence of clusters is also sequential in the sequence. An example of how to use this is
  # to output presence absence data and, cluster the resulting 1/0 matrix with a program such as jExpress and
  # take the cluster ids based on the order after clustering as input for this method. This would determine if
  # any of the blocks of clusters derived from clustering are formed of contiguous blocks of genes e.g a 
  # genomic island
  # @param [Hash] options A hash of options
  # @option options [String] :db_location The path and name of the core-access database
  # @option options [String] :output_filepath The path of the output file that wiil be created
  # @option options [Array] :cluster_order An array of cluster ids that will be retrieved in order to determine
  # the genee order in each strain relative to the cluster order specified in this option
  def output_gene_order(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'


    options = MethodArgumentParser::Parser.check_options options  do
      option :db_location, :required => true, :type => :string
      option :output_filepath, :required => true, :type => :string
    end

    output_file = File.open(options[:output_filepath], "w")

    test_connection(options[:db_location])

    strains = Strain.all
    output_file.puts "\t#{strains.map{|strain| strain.name}.join("\t")}"

    if options[:cluster_order]
      clusters = Array.new
      options[:cluster_order].each do |cluster_id|
        clusters << Cluster.find(cluster_id)
      end
    else
      clusters = get_clusters(options[:without_core_genes])
    end

    
    counter = 0
    clusters.each do |cluster|
      counter += 1
      puts "Completed output for #{counter} clusters" if counter % 100 == 0
      output_file.print "#{cluster.id}"
      ActiveRecord::Base.transaction do
        strains.each do |strain|
          gene = Gene.joins(:strain, :clusters).where("strains.id = ? AND clusters.id = ?", strain.id, cluster.id).first
          if gene.nil?
            output_file.print "\t-"
          else
            output_file.print "\t#{gene.name}"
          end
        end
      end
      output_file.puts
    end
    output_file.close
  end
  
  # A method to output a genbank format file based on the annotation in the database. The sequence file originally
  # used to make the database and associated strain name are the required inputs. The file will be created in the
  # same location as the input sequence file
  # @param [Hash] options A hash of options
  # @option options [String] :db_location The path and name of the core-access database
  # @option options [String] :strain_name The name of the strain to be annotated
  # @option options [String] :sequence_file The path to the sequence file associated with the strain and used to
  # originally create the database
  # @option options [Boolean] :merge_contigs Whether to merge the contigs to produce just a single entry in the
  # GenBank file. The position of the contigs will be recorded in the annotation if this is specified.
  def output_genbank_files_from_database(options)
    default_options = {
      :merge_contigs  => true
    }
    options.reverse_merge!(default_options)
    
    options = MethodArgumentParser::Parser.check_options options  do
      option :db_location, :required => true, :type => :string
      option :strain_name, :required => true, :type => :string
      option :sequence_file, :required => true, :type => :string
    end
    
    test_connection(options[:db_location])
    ActiveRecord::Base.transaction do
      output_file = File.open("#{options[:output_dir]}/#{File.basename(options[:sequence_file], File.extname(options[:sequence_file]))}_annotated.gbk", "w")

      sequence_objects = *rich_sequence_object_from_file(options[:sequence_file]) # * converts to array

      strain = Strain.find_by_name(options[:strain_name])
      puts "Producing annotated genbank file for #{strain.name}"
      if options[:merge_contigs]
        features = Array.new
        offset = 0
        if sequence_objects.size > 1
          sequence_objects.each do |sequence_object|
            location = "#{1 + offset}..#{sequence_object.seq.length + offset}"
            features << Bio::Feature.new('contig', location)
            offset += sequence_object.seq.length
          end
        end
        combined_sequence = ""
        sequence_objects.each do |sequence_object|
          combined_sequence += sequence_object.seq
        end
        genes = strain.genes.sort{|x,y| x.location.match(/\d+/).to_s.to_i <=> y.location.match(/\d+/).to_s.to_i}
        genes.each do |gene|
          cluster_representative_for_gene = gene.clusters.where("is_parent_cluster = ?", false).first.representative
          features << create_bio_feature(gene, cluster_representative_for_gene.annotations, options[:merge_contigs])
        end
        write_bio_sequence(:sequence => combined_sequence, :features => features, :entry_id => options[:strain_name], :definition => options[:strain_name], :output_file => output_file)
      else
        genes = strain.genes.sort{|x,y| x.name.to_i <=> y.name.to_i}
        previous_relative_location = 0
        sequence_index = 0
        features = Array.new
        genes.each do |gene|
          if gene.relative_location.match(/\d+/).to_s.to_i < previous_relative_location
            write_bio_sequence(:sequence => sequence_objects[sequence_index].seq, :features => features, :entry_id => sequence_objects[sequence_index].entry_id, :definition => sequence_objects[sequence_index].definition, :output_file => output_file)
            sequence_index += 1
            features = Array.new
          end
          cluster_representative_for_gene = gene.clusters.where("is_parent_cluster = ?", false).first.representative
          features << create_bio_feature(gene, cluster_representative_for_gene.annotations, options[:merge_contigs])
          previous_relative_location = gene.relative_location.match(/\d+/).to_s.to_i
        end
        write_bio_sequence(:sequence => sequence_objects[sequence_index].seq, :features => features, :entry_id => sequence_objects[sequence_index].entry_id, :definition => sequence_objects[sequence_index].definition, :output_file => output_file)
      end
      output_file.close
    end
  end

  private
  
  # A method to return an array of Clusters. It will fetch only clusters that are not super clusters. If the option
  # :without_core_genes is true then only "non-core" clusters will be returned
  # @param [Hash] options A hash of options
  # @option options [String] :db_location The path and name of the core-access database
  # @option options [Boolean] ::without_core_genes Whether or not to include core genes in when returning clusters
  def get_clusters(options)
    where_statement = "clusters.is_parent_cluster = ?"
    where_parameters = [false]

    if (options[:without_core_genes])
      core_cluster_ids = find_core_clusters(:db_location => options[:db_location]).map{|cluster| cluster.id}
      where_statement += " AND clusters.id NOT IN (?)"
      where_parameters << core_cluster_ids
    end

    where_array = [where_statement] + where_parameters

    return Cluster.where(where_array).all
  end
  
  # A method to create a new Bio::Feature based on a gene in the database. This will usually be a cluster
  # representative which has annotations. If there are no annotations then only the gene location will be 
  # recorded in the feature
  # @param [Gene] gene a gene from the core-access database
  # @param [Array] annotations An array of Bio::Qualifiers that will be used to annotate the feature
  # @param [Boolean] merge_contigs Whether to use the absolute or relative location. If merging contigs
  # then the absolute location will be used
  def create_bio_feature(gene, annotations,  merge_contigs)
    if merge_contigs
      cds = Bio::Feature.new('CDS', gene.location)
    else
      cds = Bio::Feature.new('CDS', gene.relative_location)
    end
    unless annotations.empty?
        annotations.each do |annotation|
        cds.append(Bio::Feature::Qualifier.new("#{annotation.qualifier}", "#{annotation.value}"))
      end
    end
    cds.append(Bio::Feature::Qualifier.new("note", "gene number: #{gene.name}"))
    return cds
  end
  
  # A method to write out a Bio::Sequence in Genbank format
  # @param [Hash] options A hash of options
  # @option options [String] :sequence The sequence as text
  # @option options [String] :entry_id The entry id (equivalent to accession)
  # @option options [String] :definition The definition line for the sequence. This should be descriptive
  # @option options [Array] :features An array of Bio::Features to be written out in GenBank format
  # @option options [String] :output_file The path to the file that will be written in GenBank format
  def write_bio_sequence(options)
    bio_sequence = Bio::Sequence.new(options[:sequence])
    bio_sequence.na
    bio_sequence.entry_id = options[:entry_id]
    bio_sequence.definition = options[:definition]
    bio_sequence.molecule_type = "DNA"
    bio_sequence.topology = "linear"
    bio_sequence.division = "PRO"
    bio_sequence.date_created = "#{Date.today.day}-#{Date::ABBR_MONTHNAMES[Date.today.month].upcase}-#{Date.today.year}"
    bio_sequence.features = options[:features]
    options[:output_file].puts bio_sequence.output(:genbank)
  end
end
