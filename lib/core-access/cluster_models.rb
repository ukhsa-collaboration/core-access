  require 'active_record'
  class Strain < ActiveRecord::Base
    has_many :genes

  end
  class Gene < ActiveRecord::Base
    belongs_to :strain
    has_many :annotations
    has_many  :cluster_memberships
    has_many :clusters, :through  => :cluster_memberships
  end
  class Cluster < ActiveRecord::Base
    has_many  :cluster_memberships
    has_many :genes, :through  => :cluster_memberships
    belongs_to :representative, :class_name => "Gene", :foreign_key => "representative_id"
    belongs_to :parent_cluster, :class_name => "Cluster", :foreign_key => "parent_id"
    has_many :child_clusters, :class_name => "Cluster", :foreign_key => "parent_id"

    def self.number_of_strains(cluster_id)
      joins(:genes => :strain).select("DISTINCT strains.id").where("clusters.id = ?", cluster_id).count
    end
  end
  class ClusterMembership < ActiveRecord::Base
    belongs_to :gene
    belongs_to :cluster
    # before_save :update_paralog_status
    # after_save :update_number_of_cluster_members
    
    def update_number_of_cluster_members
      cluster = self.cluster
      cluster.number_of_members = cluster.genes.size
      cluster.save
    end
  end
  class Annotation < ActiveRecord::Base
    belongs_to :gene
  end