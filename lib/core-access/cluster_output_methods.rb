module ClusterOutput
  # a method to return summary information about a cluster
  # @param [Cluster] cluster The cluster object (ActiveRecord) to return summary information for
  # @param [Array] attributes An array of strings listing the cluster attributes that should be returned in the summary
  # default is ["id", "cutoff", "is_parent_cluster", "number_of_members", "number_of_strains"])
  # @return [Hash] A hash of the cluster summary information. Descriptor name is the key and the info recorded in the value
  def cluster_summary(cluster, attributes = ["id", "cutoff", "is_parent_cluster", "number_of_members", "number_of_strains"])
    cluster_info = Hash.new
    attributes.each do |attribute|
      cluster_info[attribute] = cluster.send(attribute)
    end
    cluster_info["strain_names"] =  Strain.joins(:genes => :clusters).where("clusters.id = ?", cluster.id).map{|strain| strain.name}
  end

  # a method to output presence and absence data
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
    output_file.puts "\t#{strains.map{|strain| strain.name}.join("\t")}"

    clusters = get_clusters(options)

    counter = 0
    clusters.each do |cluster|
      ActiveRecord::Base.transaction do
        counter += 1
        puts "Completed output for #{counter} clusters" if counter % 100 == 0
        products = cluster.representative.annotations.select{|annotation| annotation.qualifier == "product"}
        cluster_descriptor = cluster.id.to_s
        cluster_descriptor += ": " + products.map{|product| product.value}.join(", ") unless products.empty?
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

  def output_genbank_files_from_database(options)
    default_options = {
      :merge_contigs  => true
    }
    options.reverse_merge!(default_options)
    
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
          features << create_bio_feature(cluster_representative_for_gene, options[:merge_contigs])
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
          features << create_bio_feature(cluster_representative_for_gene, options[:merge_contigs])
          previous_relative_location = gene.relative_location.match(/\d+/).to_s.to_i
        end
        write_bio_sequence(:sequence => sequence_objects[sequence_index].seq, :features => features, :entry_id => sequence_objects[sequence_index].entry_id, :definition => sequence_objects[sequence_index].definition, :output_file => output_file)
      end
      output_file.close
    end
  end

  private

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
  
  def create_bio_feature(gene, merge_contigs)
    if merge_contigs
      cds = Bio::Feature.new('CDS', gene.location)
    else
      cds = Bio::Feature.new('CDS', gene.relative_location)
    end
    unless gene.annotations.empty?
      gene.annotations.each do |annotation|
        cds.append(Bio::Feature::Qualifier.new("#{annotation.qualifier}", "#{annotation.value}"))
      end
    end
    cds.append(Bio::Feature::Qualifier.new("note", "gene number: #{gene.name}"))
    return cds
  end
  
  def write_bio_sequence(options)
    bio_sequence = Bio::Sequence.new(options[:sequence])
    bio_sequence.na
    bio_sequence.entry_id = options[:entry_id]
    bio_sequence.definition = options[:definition]
    bio_sequence.molecule_type = "DNA"
    bio_sequence.topology = "linear"
    bio_sequence.division = "PRO"
    bio_sequence.date_created = "18-NOV-2009"
    bio_sequence.features = options[:features]
    options[:output_file].puts bio_sequence.output(:genbank)
  end
end