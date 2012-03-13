class CDHitCluster
  attr_accessor :members

  def initialize
    self.members = Array.new
  end

  def add_member(member)
    self.members << member
  end
  
  def add_members(members)
    members.each do |member|
      add_member(member)
    end
  end

end

class ClusterMember
  attr_accessor :name
  attr_accessor :size
  attr_accessor :clustering_status

  def initialize(name,size,clustering_status)
    self.name = name
    self.size = size
    self.clustering_status = clustering_status
  end
end


class ClusterCollection
  attr_accessor :clusters # a hash where an arbritrary unique index is the key and the value is the cluster object
  attr_accessor :cluster_members # a hash where the member name is the key and the value is the cluster index
  
  def initialize
    self.clusters = Hash.new
    self.cluster_members = Hash.new
  end
end