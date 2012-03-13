require 'rubygems'
require 'bio'
# A method to guess a sequence format based on the filename or filepath supplied
# @param [String] filename_or_filepath The path to the sequence file whose format the method will try to guess
# @return [Symbol] the format :fasta, :genbank or :embl
def guess_sequence_format(filename_or_filepath)
  file_extension = File.extname(filename_or_filepath).downcase
  file_format = nil
  case file_extension
  when ".fa", ".fna", ".fas", ".fasta", ".fsa"
    file_format = :fasta
  when ".gbk", ".genbank", ".gb"
    file_format = :genbank
  when ".embl", ".emb"
    file_format = :embl
  end
  return file_format
end
# A method to create a bioruby Bio::Sequence object from a file
# @param [String] path_to_sequence The file path to the rich sequence file
# @param [Symbol] sequence_format The format of the sequence file. Either :genbank or :format
# @return [Bio::Sequence] The bioruby Bio::Sequence object
def biosequence_from_file(path_to_sequence, sequence_format = nil)
  biosequence = rich_sequence_object_from_file(path_to_sequence, sequence_format).to_biosequence
  return biosequence
end
# A method to create either a bioruby Bio::Genbank or Bio::EMBL object from a file
# @param [String] path_to_sequence The file path to the rich sequence file
# @param [Symbol] sequence_format The format of the sequence file. Either :genbank or :embl
# @return [Bio::Genbank,Bio::EMBL] Either a bioruby Bio::Genbank or Bio::EMBL object
def rich_sequence_object_from_file(path_to_sequence, sequence_format = nil)
  sequence_format = guess_sequence_format(path_to_sequence) if sequence_format.nil?
  rich_sequence_objects = Array.new
  case sequence_format
  when :genbank
    Bio::FlatFile.open(Bio::GenBank, path_to_sequence).each_entry do |entry|
      rich_sequence_objects << entry unless entry.entry_id.nil? || entry.entry_id.empty?
    end
  when :embl
    Bio::FlatFile.open(Bio::EMBL, path_to_sequence).each_entry do |entry|
      rich_sequence_objects << entry unless entry.entry_id.nil? || entry.entry_id.empty?
    end
  end
  if rich_sequence_objects.size > 1
    return rich_sequence_objects
  else
    return rich_sequence_objects.first
  end
end
# A method to create a text file representing sequence from a bioruby Bio::Sequence object
# @param [String] path_to_sequence The file path to the rich sequence file that will be created
# @param [Symbol] sequence_format The format of the sequence file to be created. Either :fasta, :genbank or :format
# @param [Bio::Sequence] biosequence The biosequence to be converted to a text file
def file_from_biosequence(path_to_sequence, sequence_format, biosequence)
  ouput_file = File.open(path_to_sequence, "w")
  ouput_file.puts biosequence.output(sequence_format)
  ouput_file.close
end
# A method to convert a sequence file from one rich format to another
# @param [String] path_to_sequence The file path to the rich sequence file that will be created
# @param [Symbol] input_sequence_format The format of the sequence file. Either :genbank or :embl
# @param [Symbol] output_sequence_format The format of the sequence file to be created. Either :fasta, :genbank or :embl
def change_sequence_format(path_to_sequence, input_sequence_format, output_sequence_format)
  biosequence = biosequence_from_file(path_to_sequence, input_sequence_format)
  case output_sequence_format
  when :genbank
    output_suffix = "gbk"
  when :embl
    output_suffix = "embl"
  when :fasta
    output_suffix = "fas"
  end
  output_filepath = "#{File.dirname(path_to_sequence)}/#{file_prefix_from_filepath(path_to_sequence)}.#{output_suffix}"
  file_from_biosequence(output_filepath, output_sequence_format, biosequence)
end
# a method to create a biosequence based on a feature extracted from a rich sequence object
# @param [Bio::Feature] cds The feature which will be used to generate the biosequence
# @param [Bio::(Genbank::EMBL)] rich_sequence_object The rich sequence object from which the feature 
# was derived
def create_biosequence_from_feature(feature,rich_sequence_object, protein = false)
  location_string = location_string_from_feature(feature)
  cds_sequence = rich_sequence_object.seq.splice(location_string)
  cds_sequence = cds_sequence.translate(1,11) if protein
  query_sequence = Bio::Sequence.new(cds_sequence.to_s)
end
# a method to extract the location string that can be used to splice out sequence or for other purposes
# @param [Bio::Feature] cds The feature which will be used to generate the biosequence
# @return [String] The location string of the feature
def location_string_from_feature(feature)
  location =  feature.locations.first
  location_string = location.range.to_s
  if location.strand == -1 #negative strand
    location_string = "complement(#{location_string})"
  end
  location_string
end
# a method to convert a sequence a rich format to fasta
# @param [String] path_to_sequence path to the input sequence in genbank or embl format
# @param [String] output_sequence path to the output sequence or a file object to write to. Could be a TempFile object
def rich_sequence_to_fasta(path_to_sequence, output_sequence)
  require 'tempfile'
  sequence_format = guess_sequence_format(path_to_sequence)
  case sequence_format
  when :genbank
    rich_sequence_object = Bio::FlatFile.open(Bio::GenBank, path_to_sequence).first
  when :embl
    rich_sequence_object = Bio::FlatFile.open(Bio::EMBL, path_to_sequence).first
  end
  biosequence = rich_sequence_object.to_biosequence
  
  case output_sequence
  when String
    file_from_biosequence(output_sequence, :fasta, biosequence)
  when Tempfile
    file_from_biosequence(output_sequence.path, :fasta, biosequence)
  end
end

# a method to reverse complement a rich sequence
# @param [String] path_to_sequence Path to the rich sequence file
# @param [Symbol] sequence_format The format of the rich sequence either :embl or :genbank
# @param [Symbol] output_sequence_format The format inwhich to output the  sequence either :embl or :genbank, by default will be the same as input
def reverse_complement_rich_sequence(path_to_sequence, sequence_format = nil, output_sequence_format = nil)
  biosequence = biosequence_from_file(path_to_sequence, sequence_format)
  biosequence.reverse_complement!

  # now to reverse complement the features
  reverse_complemented_features = Array.new
  sequence_length = biosequence.to_s.length
  source_feature = nil

  biosequence.features.each do |feature|
    location = feature.locations.locations.first
    if feature.feature == "source"
      source_feature = feature
    else
      new_start = sequence_length - location.to + 1
      new_end =  sequence_length - location.from + 1
      location.from = new_start
      location.to = new_end
      location.strand *= -1 # reverse the strand

      feature_locations = feature.locations
      feature_locations.locations = [location]
      feature.position = feature_locations.to_s # the location of a feature is based on it's position string, not a Bio:Locations object which is created on the fly (why I have no idea!)
      reverse_complemented_features.unshift(feature)
    end
  end
  reverse_complemented_features.unshift(source_feature)
  biosequence.features = reverse_complemented_features

  output_sequence_format ||= guess_sequence_format(path_to_sequence)
  case output_sequence_format
  when :genbank
    output_suffix = "gbk"
  when :embl
    output_suffix = "embl"
  end
  reverse_complemented_file = File.open("#{File.dirname(path_to_sequence)}/#{file_prefix_from_filepath(path_to_sequence)}_rc.#{output_suffix}", "w")
  
  reverse_complemented_file.puts biosequence.output(output_sequence_format)
  reverse_complemented_file.close
end

# a method to change the origin of a sequence
# @param [String] path_to_sequence Path to the rich sequence file
# @param [Integer] new_origin The sequence position that is the new origin of the sequence
# @param [Symbol] sequence_format The format of the rich sequence either :embl or :genbank
# @param [Symbol] output_sequence_format The format inwhich to output the  sequence either :embl or :genbank, by default will be the same as input
def change_origin_of_rich_sequence(path_to_sequence, new_origin, sequence_format = nil, output_sequence_format = nil)
  biosequence = biosequence_from_file(path_to_sequence, sequence_format)
  sequence_length = biosequence.to_s.length
  # shift the sequence
  first_sequence_part = biosequence.splice((1..new_origin-1).to_s)
  last_sequence_part = biosequence.splice((new_origin..sequence_length).to_s)
  biosequence.seq = last_sequence_part + first_sequence_part

  # now to shift those features
  five_prime_shifted_features = Array.new
  three_prime_shifted_features = Array.new
  source_feature = nil

  biosequence.features.each do |feature|
    location = feature.locations.locations.first
    if feature.feature == "source"
      source_feature = feature
    else
      if location.from >= new_origin
        offset = -new_origin + 1
      else
        offset = sequence_length - new_origin + 1
      end
      new_start = location.from + offset
      new_end =  location.to + offset
      location.from = new_start
      location.to = new_end

      feature_locations = feature.locations
      feature_locations.locations = [location]
      feature.position = feature_locations.to_s # the location of a feature is based on it's position string, not a Bio:Locations object which is created on the fly (why I have no idea!)
      if offset < 0 # negative offset means that this feature needs to go into new 5_prime_shifted_features
        five_prime_shifted_features << feature
      else
        three_prime_shifted_features << feature
      end
    end
  end
  five_prime_shifted_features.unshift(source_feature)
  biosequence.features = five_prime_shifted_features + three_prime_shifted_features

  output_sequence_format ||= guess_sequence_format(path_to_sequence)
  case output_sequence_format
  when :genbank
    output_suffix = "gbk"
  when :embl
    output_suffix = "embl"
  end
  reoriented_file = File.open("#{File.dirname(path_to_sequence)}/#{file_prefix_from_filepath(path_to_sequence)}_reoriented.#{output_suffix}", "w")
  
  puts  biosequence.output(output_sequence_format)
  # reoriented_file.puts biosequence.output(output_sequence_format)
  reoriented_file.close
end
# A method to find the feature at a given position
def find_features_at_position(biosequence, *positions)
  positions.flatten!
  features = Array.new
  biosequence.features.each do |feature|
    out_of_range = 0
    positions.each do |position|
      if position >= feature.locations.span.first && position <= feature.locations.span.last
        features << feature unless feature.feature == "source" || features.include?(feature)
      elsif feature.locations.span.first > position
        out_of_range += 1
      end
    end
    break if out_of_range == positions.size
  end
  return features
end

# # # A method to get the file prefix based on a file path
# @param [String] filepath The filepath
# @return [String] The file prefix
def file_prefix_from_filepath(filepath)
  File.basename(filepath).sub(/#{File.extname(filepath)}$/,'')
end
  