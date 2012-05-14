require 'rubygems'
require 'core-access/cluster_create_methods'
extend ClusterCreate
require 'core-access/cluster_annotate_methods'
extend ClusterAnnotate
require 'core-access/cluster_search_methods'
extend ClusterSearch
require 'core-access/cluster_output_methods'
extend ClusterOutput
require 'core-access/cluster_database'
extend ClusterDB
require 'core-access/cluster_models'

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
  if choices[:sequence_file_list]
    file_info = File.read(choices[:sequence_file_list]).split("\n")
  elsif choices[:genbank_output_strain_list_file]
    file_info = File.read(choices[:genbank_output_strain_list_file]).split("\n")
  end
  sequence_files = file_info.map{|line| line.split("\t").first}
  # optionally append a common sequence file directory
  sequence_files.map!{|sf| "#{choices[:sequence_file_dir]}/#{sf}"} if choices[:sequence_file_dir]

  strain_names = file_info.map{|line| line.split("\t").last}
  return sequence_files, strain_names
end