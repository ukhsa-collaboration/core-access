module ClusterOutput

  def cluster_summary(cluster, attributes = ["id", "cutoff", "is_parent_cluster", "number_of_members", "number_of_strains"])
    cluster_info = Array.new
    attributes.each do |attribute|
      cluster_info << cluster.send(attribute)
    end
    cluster_info <<  Strain.joins(:genes => :clusters).where("clusters.id = ?", cluster.id).map{|strain| strain.name}
  end


  def output_gene_presence_absence(options)
    require 'core-access/cluster_database'
    extend ClusterDB
    require 'core-access/cluster_models'

    default_options = {
      :without_core_genes  => false,
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

    clusters = get_clusters(options[:without_core_genes])

    counter = 0
    clusters.each do |cluster|
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
    test_connection(options[:db_location])

    output_file = File.open("#{options[:output_dir]}/#{File.basename(options[:sequence_file], File.extname(options[:sequence_file]))}_annotated.gbk", "w")
    combined_sequence = ""
    features = Array.new

    sequence_objects = *rich_sequence_object_from_file(options[:sequence_file]) # * converts to array
    offset = 0
    if sequence_objects.size > 1
      sequence_objects.each do |sequence_object|
        location = "#{1 + offset}..#{sequence_object.seq.length + offset}"
        features << Bio::Feature.new('contig', location)
        offset += sequence_object.seq.length
      end
    end
    sequence_objects.each do |sequence_object|
      combined_sequence += sequence_object.seq
    end

    strain = Strain.find_by_name(options[:strain_name])
    genes = strain.genes.sort{|x,y| x.location.match(/\d+/).to_s.to_i <=> 
      y.location.match(/\d+/).to_s.to_i}
    genes.each do |gene|
      cds = Bio::Feature.new('CDS', gene.location)
      unless gene.annotations.empty?
        gene.annotations.each do |annotation|
          cds.append(Bio::Feature::Qualifier.new("#{annotation.qualifier}", "#{annotation.value}"))
        end
      end
      cds.append(Bio::Feature::Qualifier.new("note", "gene number: #{gene.name}"))
      features << cds
    end
    bio_sequence = Bio::Sequence.new(combined_sequence)
    bio_sequence.na
    bio_sequence.entry_id = options[:strain_name]
    bio_sequence.definition = options[:strain_name]
    bio_sequence.molecule_type = "DNA"
    bio_sequence.topology = "linear"
    bio_sequence.division = "PRO"
    bio_sequence.date_created = "18-NOV-2009"
    bio_sequence.features = features
    output_file.puts bio_sequence.output(:genbank)
    output_file.close
  end

  private

  def get_clusters(options, without_core_genes = false)
    where_statement = "clusters.is_parent_cluster = ?"
    where_parameters = [false]

    if (without_core_genes)
      core_cluster_ids = find_core_clusters(:db_location => options[:db_location]).map{|cluster| cluster.id}
      where_statement += " AND clusters.id NOT IN (?)"
      where_parameters << core_cluster_ids
    end

    where_array = [where_statement] + where_parameters

    return Cluster.where(where_array).all
  end
end