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

def print_help_for(command)
  choice_options = Choice.class_variable_get("@@options").dup
  choice_option = nil
  common_options = Array.new

  until choice_option.class == String && choice_option =~ /Arguments common to all commands/ do
    choice_option = choice_options.shift
  end

  until (choice_option.class == String && choice_option =~ /Arguments for/) || choice_option.nil? do
    choice_option = choice_options.shift
    common_options << choice_option.last if choice_option.class == Array
  end

  until choice_option.class == String && choice_option =~ /Arguments for #{command}/ do
    choice_option = choice_options.shift
  end

  command_options = Hash.new
  if command.empty?
    commands = ["create", "annotate", "search", "output"]
  else
    commands = [command]
  end

  commands.each do |command|
    command_options[command] = Array.new
    until (choice_option.class == String && choice_option !~ /Arguments for #{command}/ && choice_option =~ /Arguments for/) || choice_option.nil? do
      choice_option = choice_options.shift
      command_options[command] << choice_option.last if choice_option.class == Array
    end
  end

  puts "Usage: core-access #{command} OPTIONS"
  puts "Arguments common to all commands:"
  common_options.each do |common_option|
    print_option(common_option)
  end
  puts
  commands.each do |command|
    puts "Options for the core-access #{command} command:"
    command_options[command].each do |command_option|
      print_option(command_option)
    end
    puts
  end
end

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

private
def print_option(option)
  # Make this easier on us
  short = option.short
  long = option.long
  line = ''

  # Get the short part.
  line << sprintf("%6s", short)
  line << sprintf("%-2s", (',' if short && long))

  # Get the long part.
  line << sprintf("%-29s", long)

  # Print what we have so far
  print line

  # If there's a desc, print it.
  if option.desc
    # If the line is too long, spill over to the next line
    if line.length > 37
      puts           
      print " " * 37
    end

    puts option.desc.shift
    
    # If there is more than one desc line, print each one in succession
    # as separate lines.
    option.desc.each do |desc| 
      puts ' '*37 + desc
    end

  else
    # No desc, just print a newline.
    puts 
  end
end