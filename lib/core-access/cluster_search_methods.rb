module ClusterSearch
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

    sql_statement += " AS cl" if options[:unique_to_subset]

    sql_statement += " WHERE number_of_strains = #{strain_names.size}" if options[:unique]
    if options[:unique_to_subset]
      sql_statement += " GROUP BY id HAVING COUNT(*) = cl.number_of_strains"
    else
      sql_statement += " GROUP BY id HAVING COUNT(*) = #{strain_names.size}"
    end
    return Cluster.find_by_sql(sql_statement)
  end

  def find_unique_clusters(options) # this method finds genes that are shared by all strains supplied and found nowhere else
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
end