@@test_dir = File.dirname(__FILE__)
require 'helper'
require 'stringio'
require 'fileutils'
require "#{@@test_dir}/../lib/core-access" # require lib/core-access so changes are picked up without rake install
include Clustering

class TestCoreAccess < Test::Unit::TestCase
  @@clustering_cutoffs = [99,98,95,90,85]
  # context "core-access test" do

  #   should "find the sequence file list file" do
  #     assert File.exists?("#{@@test_dir}/fasta_sequence_file_list.txt")
  #   end
  #   should "find all files in sequence list" do
  #     file_exists_array = Array.new
  #     File.read("#{@@test_dir}/fasta_sequence_file_list.txt").split("\n").map{|line| line.split("\t").first}.each do |file|
  #       file_exists_array << File.exists?("#{@@test_dir}/#{file}")
  #     end
  #     assert_equal(false, file_exists_array.include?(false))
  #   end
  #   should "find the reference file" do
  #     assert File.exists?("#{@@test_dir}/test_data/ref.gbk")
  #   end
  # end

  context "core-access methods" do
  #   should "find the correct strain names when supplied" do
  #     sequence_files, strain_names = extract_file_and_strain_names_from_file_list(:sequence_file_list => "#{@@test_dir}/fasta_sequence_file_list.txt")
  #     assert_equal ["strain1", "strain2"], strain_names
  #   end
    # should "create a multi fasta file from the input files" do
    #   FileUtils.rm("/tmp/cds_proteins.fas", :force => true)
    #   sequence_files, strain_names = extract_file_and_strain_names_from_file_list(:sequence_file_list => "#{@@test_dir}/fasta_sequence_file_list.txt")
    #   sequence_files.map!{|sf| "#{@@test_dir}/#{sf}"}
    #   create_cds_multi_fasta_file(:training_sequence_path => "#{@@test_dir}/test_datas/ref.gbk", :root_folder  => "/tmp", :cds_multi_fasta_file => "/tmp/cds_proteins.fas", :sequence_files => sequence_files)
    #   assert File.exists?("/tmp/cds_proteins.fas")
    # end
    # should "create cd-hit output from a multifasta file" do
    #   FileUtils.rm("/tmp/cdhit_clusters*", :force => true)
    #   make_heirachical_clusters(:root_folder  => "/tmp", :cds_multi_fasta_file => "#{@@test_dir}/test_data/cds_proteins.fas", :clustering_cutoffs => @@clustering_cutoffs, :cdhit_dir => "/usr/local/cdhit")
    #   assert File.exists?("/tmp/cdhit_clusters-85")
    # end
    # should "create a cluster database from cd-hit output" do
    #   FileUtils.rm("/tmp/test.sqlite", :force => true)
    #   sequence_files, strain_names = extract_file_and_strain_names_from_file_list(:sequence_file_list => "#{@@test_dir}/fasta_sequence_file_list.txt")
    #   @@clustering_cutoffs.each do |clustering_cutoff|
    #     FileUtils.cp("#{@@test_dir}/test_data/cdhit_clusters-#{clustering_cutoff}.clstr", "/tmp/")
    #   end
    #   # create database
    #   make_db_clusters(:root_folder => "/tmp", :db_location => "/tmp/test.sqlite", :clustering_cutoffs => @@clustering_cutoffs, :strain_names => strain_names)
    #   assert_equal 2, Cluster.where(:id => 1).first.genes.size
    #   assert_equal 1, Cluster.where(:id => Cluster.count).first.genes.size
    # end

    # should "create representatives in clusters" do
    #   FileUtils.rm("/tmp/test.sqlite", :force => true)
    #   FileUtils.cp("#{@@test_dir}/test_data/test_without_representatives.sqlite", "/tmp/test.sqlite")
    #   # add represenatives to clusters
    #   sequence_files, strain_names = extract_file_and_strain_names_from_file_list(:sequence_file_list => "#{@@test_dir}/genbank_sequence_file_list.txt")
    #   sequence_files.map!{|sf| "#{@@test_dir}/#{sf}"}
    #   add_representative_sequences_to_cluster(:root_folder => "/tmp", :db_location => "/tmp/test.sqlite", :sequence_files => sequence_files, :strain_names => strain_names)
    #   assert_equal @@representative_of_cluster1_sequence, Cluster.where(:id => 1).first.representative.sequence
    #   assert_equal @@representative_of_cluster2_sequence, Cluster.where(:id => Cluster.count).first.representative.sequence
    # end

    # should "make super clusters" do
    #   FileUtils.rm("/tmp/test.sqlite", :force => true)
    #   FileUtils.cp("#{@@test_dir}/test_data/test_with_representatives.sqlite", "/tmp/test.sqlite")
    #   # add represenatives to clusters
    #   make_super_clusters(:cdhit_dir => "/usr/local/cdhit", :cutoff => 65, :db_location => "/tmp/test.sqlite")
    #   puts Cluster.count
    #   assert_equal true, Cluster.where(:id => Cluster.count).first.is_parent_cluster
    #   assert_equal 65, Cluster.where(:id => Cluster.count).first.cutoff
    #   assert_equal 2, Cluster.where(:id => Cluster.count).first.number_of_members
    #   assert_equal 1, Cluster.where(:id => Cluster.count).first.number_of_strains
    # end

    # should "find correct number of core genes" do
    #   FileUtils.rm("/tmp/test.sqlite", :force => true)
    #   FileUtils.cp("#{@@test_dir}/test_data/test_with_superclusters.sqlite", "/tmp/test.sqlite")
    #   # find core genes
    #   core_clusters = find_core_clusters(:db_location => "/tmp/test.sqlite")
    #   assert_equal 326, core_clusters.size
    # end

    should "annotate clusters" do
      FileUtils.rm("/tmp/test.sqlite", :force => true)
      FileUtils.cp("#{@@test_dir}/test_data/test_with_superclusters.sqlite", "/tmp/test.sqlite")
      annotate_clusters(:db_location => "/tmp/test.sqlite", :root_folder => "/tmp", :reference_genomes => ["#{@@test_dir}/test_data/ref_seq1.gbk","#{@@test_dir}/test_data/ref_seq2.gbk" ], :microbial_genomes_blast_db => "/Volumes/DataRAID/blast_databases/microbial_genomes", :blast_dir => "/usr/local/blast/bin")
    end


  end

  context "core-access" do
  #     should "fail to find glimmer3 in /tmp" do
  #     output = `#{@@test_dir}/../bin/core-access create -db test.sqlite3 -od /tmp -fl #{@@test_dir}/fasta_sequence_file_list.txt -fd #{@@test_dir} -gd /tmp 2>&1`
  #     assert output =~ /Can not find glimmer3/
  #   end
  #   should "fail to find cd-hit on /tmp" do
  #     output = `#{@@test_dir}/../bin/core-access create -db test.sqlite3 -od /tmp -fl #{@@test_dir}/genbank_sequence_file_list.txt -gd /usr/local/glimmer -fd #{@@test_dir} -cd /tmp 2>&1`
  #     assert output =~ /Can not find cd-hit/
  #   end
  #   should "create a multiple fasta file, run cd-hit and make a database when specifying the create command and correct options" do
  #     FileUtils.rm("/tmp/cds_proteins.fas", :force => true)
  #     FileUtils.rm("/tmp/cdhit_clusters*", :force => true)
  #     FileUtils.rm("/tmp/test.sqlite", :force => true)
  #     output = `#{@@test_dir}/../bin/core-access create -db test.sqlite -fl #{@@test_dir}/fasta_sequence_file_list.txt -fd #{@@test_dir} -gd /usr/local/glimmer -cd /usr/local/cdhit -sc 65 -od /tmp 2>&1`
  #     puts output
  #     # check that the multi fasta file creation worked
  #     assert File.exists?("/tmp/cds_proteins.fas")
  #     # check that cd-hit ran OK
  #     assert File.exists?("/tmp/cdhit_clusters-85")
  #     # check that database was created with correct values
  #     assert File.exists?("/tmp/test.sqlite")
  #     require 'core-access/cluster_database'
  #     extend ClusterDB
  #     require 'core-access/cluster_models'
  #     make_db_connection("/tmp/test.sqlite")
  #     assert_equal 2, ::Cluster.where(:id => 1).first.genes.size
  #     assert_equal 1, ::Cluster.where(:id => Cluster.where(:is_parent_cluster => false).count).first.genes.size
  #     # check representatives are in the db
  #     assert_equal @@representative_of_cluster1_sequence, Cluster.where(:id => 1).first.representative.sequence
  #     assert_equal @@representative_of_cluster2_sequence, Cluster.where(:id => Cluster.where(:is_parent_cluster => false).count).first.representative.sequence
  #     # check super clusters are in database
  #     assert_equal true, Cluster.where(:id => Cluster.count).first.is_parent_cluster
  #     assert_equal 65, Cluster.where(:id => Cluster.count).first.cutoff
  #     assert_equal 3, Cluster.where(:id => Cluster.count).first.number_of_members
  #     assert_equal 2, Cluster.where(:id => Cluster.count).first.number_of_strains
  #   end
    # should "find correct number of core genes" do
    #   FileUtils.rm("/tmp/test.sqlite", :force => true)
    #   FileUtils.cp("#{@@test_dir}/test_data/test_with_superclusters.sqlite", "/tmp/test.sqlite")
    #   output = `#{@@test_dir}/../bin/core-access search -cg -db /tmp/test.sqlite`
    #   puts output
    # end
  end
end
