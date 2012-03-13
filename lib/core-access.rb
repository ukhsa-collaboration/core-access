require 'rubygems'
require 'core-access/clustering_methods'


def determine_glimmer_directory(choices)
  if choices[:glimmer_dir]
    glimmer_dir = choices[:glimmer_dir]
  else
    glimmer_dir = which("glimmer3")
    glimmer_dir = File.dirname(glimmer_dir).sub(/\/bin/,'') unless glimmer_dir.nil? 
  end
  return glimmer_dir
end

def extract_file_and_strain_names_from_file_list(choices)
  file_info = File.read(choices[:sequence_file_list]).split("\n")
  sequence_files = file_info.map{|line| line.split("\t").first}
  sequence_files.map!{|sf| "#{choices[:sequence_file_dir]}/#{sf}"} if choices[:sequence_file_dir]
  strain_names = file_info.map{|line| line.split("\t").last}
  return sequence_files, strain_names
end