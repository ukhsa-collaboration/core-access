module ClusterDB

  def make_db_connection(db_location)
    require 'active_record'
      ActiveRecord::Base.establish_connection(
        :adapter => "sqlite3",
        :database => db_location
      )
  end

  def test_connection(db_location = nil)
    require 'active_record'
    begin
      ActiveRecord::Base.connection
    rescue ActiveRecord::ConnectionNotEstablished
      unless db_location.nil?
        make_db_connection(db_location)
      else
        raise "You need to establish a connection to the database using the make_db_connection(db_location) method"
      end
    end
  end

  def make_cluster_db_schema
    require 'active_record'
    ActiveRecord::Schema.define do
      unless table_exists? :strains
        create_table :strains do |t|
          t.column :name, :string
          t.column :description, :string
        end
      end

      unless table_exists? :genes
        create_table :genes do |t|
          t.column :name,      :string
          t.column :aa_length, :string
          t.column :location,  :string
          t.column :relative_location, :string # relative location is the location within the original contig
          t.column :sequence, :string
          t.column :strain_id, :integer
        end
      end

      unless table_exists? :clusters
        create_table :clusters do |t|
          t.column :name, :string
          t.column :parent_id, :integer
          t.column :is_parent_cluster, :boolean
          t.column :cutoff, :integer
          t.column :representative_id, :integer
          t.column :number_of_members, :integer
          t.column :number_of_strains, :integer
          t.column :contains_paralogs, :boolean
        end
      end

      unless table_exists? :cluster_memberships
        create_table :cluster_memberships do |t|
          t.column :cluster_id, :integer
          t.column :gene_id, :integer
          t.column :status, :string
        end
      end
      
      unless table_exists? :annotations
        create_table :annotations do |t|
          t.column  :qualifier, :string
          t.column  :value, :string
          t.column  :gene_id, :integer
        end
      end
       
      # indices
      unless index_exists? :genes, :name
        add_index :genes, :name
      end
      unless index_exists? :genes, :strain_id
        add_index :genes, :strain_id
      end
      unless index_exists? :clusters, :number_of_members
        add_index :clusters, :number_of_members
      end
      unless index_exists? :clusters, :number_of_strains
        add_index :clusters, :number_of_strains
      end
      unless index_exists? :clusters, :contains_paralogs
        add_index :clusters, :contains_paralogs
      end
      unless index_exists? :cluster_memberships, :cluster_id
        add_index :cluster_memberships, :cluster_id
      end
      unless index_exists? :cluster_memberships, :gene_id
        add_index :cluster_memberships, :gene_id
      end
      unless index_exists? :annotations, :gene_id
        add_index :annotations, :gene_id
      end
    end
  end
end
