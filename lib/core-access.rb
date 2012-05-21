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
    glimmer_dir = which("glimmer3") # method found in utils/file_utils
    glimmer_dir = File.dirname(glimmer_dir).sub(/\/bin/,'') unless glimmer_dir.nil? 
  end
  return glimmer_dir
end

def determine_cdhit_directory(choices)
  if Choice.choices[:cdhit_dir]
    cdhit_dir = Choice.choices[:cdhit_dir]
  else
    cdhit_dir = which("cdhit")
    cdhit_dir = File.dirname(cdhit_dir) unless cdhit_dir.nil? 
  end
  return cdhit_dir
end

def extract_file_and_strain_names_from_file_list(choices)
  if choices[:sequence_file_list]
    file_info = File.read(choices[:sequence_file_list]).split("\n")
  elsif choices[:genbank_strain_file_list]
    file_info = File.read(choices[:genbank_strain_file_list]).split("\n")
  end
  sequence_files = file_info.map{|line| line.split("\t").first}
  # optionally append a common sequence file directory
  sequence_files.map!{|sf| "#{choices[:sequence_file_dir]}/#{sf}"} if choices[:sequence_file_dir]

  strain_names = file_info.map{|line| line.split("\t").last}
  return sequence_files, strain_names
end

def silence &block
  default_stdout = $stdout
  $stdout = StringIO.new
  yield
  $stdout.string
ensure
  $stdout = default_stdout
end

def processing_indicator(indicator_pos)
  # move cursor to beginning of line
  cr = "\r"           
  # ANSI escape code to clear line from cursor to end of line
  # "\e" is an alternative to "\033"
  # cf. http://en.wikipedia.org/wiki/ANSI_escape_code
  clear = "\e[0K"     

  # reset lines
  reset = cr + clear
  indicator_characters = [ "|", "/", "-", "\\" ]

  # 7 turns on reverse video mode, 31 red , 32 green ...
  n = 32

  case indicator_pos
  when 0..3
    prefix = "#{reset}\e[#{n};1m"   
    print "#{prefix}#{indicator_characters[indicator_pos]}"
  else
    print "\e[0m"
    print "#{reset}"
  end
  $stdout.flush
end