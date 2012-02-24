##
# Two additional methods to allow for reverse merging two hashes where the keys in the calling hash take precedence over those in the other_hash.
# This is particularly useful for initializing an option hash with default values:
# 
# @author Anthony Underwood
class Hash
  # Performs the opposite of merge, with the keys and values from the first hash taking precedence over the second.
  # @param [Hash] other_hash A second hash which will be merged into the first but with he first hash taking precedence over this one
  def reverse_merge(other_hash)
    other_hash.merge(self)
  end
  
  # Performs the opposite of merge, with the keys and values from the first hash taking precedence over the second. Modifies the receiver in place.
  # @param [Hash] other_hash A second hash which will be merged into the first but with he first hash taking precedence over this one
  def reverse_merge!(other_hash)
    replace(reverse_merge(other_hash))
  end
end