module ClusterOutput

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
      strain.each do |strain|
        gene = Gene.joins(:strain, :clusters).where("strains.id = ? AND clusters.id = ?", strain.id, cluster.id).first
        if gene.nil?
          output_file.print "\t-"
        else
          output_file.print "\t#{gene.name}"
        end
      end
    end
  end

  private

  def get_clusters(without_core_genes = false)
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