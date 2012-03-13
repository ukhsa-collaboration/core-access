##
# A set of methods to perform reciprocal blast hits either against a local or remote blast database
# 
# @author Anthony Underwood
module ReciprocalHit
  require 'utils/method_argument_parser'
  require 'utils/hash_reverse_merge'
  require 'json'
  
  require 'bioutils/blast'

  include MethodArgumentParser
  include Blast

  # A method to perform a reciprocal blast search against a local blast database
  # @param [Hash] options The options for the method
  # @option options [Bio::Sequence] :query_sequence A Bio::Sequence where the entry_id has been set to the id of the sequence.
  # This should be the same as has been set in the local blast database containing the query sequence
  # @option options [String] :query_sequence_blast_database The path and name of the local blast database that contains the query sequence so that 
  # the 'reverse' blast can hopfully find the query sequence
  # @option options [String] :subject_blast_database The path and name of the local blast database from which you are trying to find a reciprocal hit
  # @option options [String] :path_to_blast_executable The path to the blast executable, default is "/usr/local/blast/bin/blastall",
  # @option options [String] :blast_options Blast options e.g default is '-e 1e-10' (see http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastall/blastall_node23.html)
  #  @option options [Boolean] :self_search search self for reciprocal hit (paralogs)
  #  @option options [String] :blast_program Which blast program to use, e.g blastn, blastp
  # @return [Boolean, Array] This method will return true or false depending on whether a reciprocal blast hit could be found in the database.
  # If a hit was found then the accession/id and descripton wil be returned also. If not then the reason for failure will be returned
  def find_local_reciprocal_hit(options)
    default_options = {
      :path_to_blast_executable => "/usr/local/blast/bin/blastall",
      :blast_options => '-e 1e-20 -F F -b 10 -v 10',
      :fastacmd_dir => "/usr/local/blast/bin/",
      :self_search => false,
      :blast_program => 'blastn',
      :percent_identity_cutoff => 90,
      :minimum_hit_length => 85,
      :display_metrics => false
    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :query_sequence, :required => true, :type => 'Bio::Sequence'
      option :query_sequence_blast_database,  :required => true, :type => :string
      option :subject_blast_database, :required => true, :type => :string
    end
    
    report = blast_sequence_locally_without_bioruby(:sequence => options[:query_sequence].seq, :path_to_blast_executable => options[:path_to_blast_executable], :blast_program => options[:blast_program], :blast_database => options[:subject_blast_database], :blast_options => options[:blast_options])
    if report.nil?
      first_subject_hit = nil
    elsif options[:self_search]
      first_subject_hit = report.hits[1] # second hit because first will be itself
    else
      first_subject_hit = report.hits.first
    end
    # puts first_subject_hit.definition
    if first_subject_hit.nil?
      return false, "No reciprocal hit by local blast: no hit in subject database"
    else
      if significant_hit?(:hit => first_subject_hit, :query_sequence => options[:query_sequence], :percent_identity_cutoff => options[:percent_identity_cutoff], :minimum_hit_length => options[:minimum_hit_length], :display_metrics => options[:display_metrics])
        if options[:blast_program] == "blastn"
          protein_option = "F"
        else
          protein_option = "T"
        end
        hit_sequence_details = retrieve_sequence(:blast_database => options[:subject_blast_database], :search_string => first_subject_hit.accession, :fastacmd_dir => options[:fastacmd_dir], :protein_option => protein_option)
        hit_sequence_object = Bio::FastaFormat.new(hit_sequence_details)
        # do reciprocal blast
        reciprocal_hit_found, number_of_reciprocal_hits_found, message = perform_reciprocal_blast(:query_sequence_for_reciprocal_blast => hit_sequence_object, :original_query_sequence_id => options[:query_sequence].entry_id, :path_to_blast_executable => options[:path_to_blast_executable],:blast_program => options[:blast_program], :query_sequence_blast_database => options[:query_sequence_blast_database], :percent_identity_cutoff => options[:percent_identity_cutoff], :minimum_hit_length => options[:minimum_hit_length], :blast_options => options[:blast_options])
        if reciprocal_hit_found
          # puts "number_of_reciprocal_hits_found: #{number_of_reciprocal_hits_found}"
          if number_of_reciprocal_hits_found == 1 # just one reciprocal hit to original
            return true, first_subject_hit.hit_id, JSON.parse(first_subject_hit.definition)
          else # multiple significant reciprocal hits
            return false, first_subject_hit.hit_id, JSON.parse(first_subject_hit.definition)
          end
        else
          return false, message
        end
      else
        return false, "No reciprocal hit by local blast: the hit in the subject database was not significant" # no sig hit so not a recipiprocal hit
      end
    end
  end
  # A method to perform a reciprocal blast search against a remote blast database
  # @param [Hash] options The options for the method
  # @option options [Bio::Sequence] :query_sequence A Bio::Sequence where the entry_id has been set to the id of the sequence.
  # This should be the same as has been set in the local blast database containing the query sequence
  # @option options [String] :query_sequence_blast_database The path and name of the local blast database that contains the query sequence so that 
  # the 'reverse' blast can hopfully find the query sequence
  # @option options [String] :blast_program The blast program (blastn is the default)
  # @option options [String] :blast_database The remote database to be searched ('nr' by default)
  # @option options [Hash] :blast_options A hash of options to fine tune the blast search e.g :expect => '1e-3', :filter => 'T',
  # @option options [Hash] :embl_objects An optional hash of embl objects so that it is not necessary
  # to download an embl object every time if it has aleady been downloaded
  # @option options [String] :path_to_blast_executable The path to the blast executable, default is "/usr/local/blast/bin/blastall", this is ncessary for the reverse blast against the query sequence database
  # @return [Boolean, Array] This method will return true or false depending on whether a reciprocal blast hit could be found in the database.
  # If a hit was found then the accession/id and descripton wil be returned also. If not then the reason for failure will be returned
  def find_remote_reciprocal_hit(options)
    default_options = {
      :blast_program => "blastn",
      :blast_database => "nr",
      :blast_options => "",
      :path_to_blast_executable => "/usr/local/blast/bin/blastall",
      :embl_objects => {},
      :percent_identity_cutoff => 90,
      :minimum_hit_length => 85

    }
    options.reverse_merge!(default_options)

    options = MethodArgumentParser::Parser.check_options options  do
      option :query_sequence, :required => true, :type => 'Bio::Sequence'
      option :query_sequence_blast_database,  :required => true, :type => :string
    end
    report = nil
    blast_attempts = 0
    while report.nil? || blast_attempts < 6
      blast_attempts += 1
      report = blast_sequence_remotely(:sequence => options[:query_sequence].seq, :blast_program => options[:blast_program], :blast_database => options[:blast_database], :blast_options => options[:blast_options])
    end
    if report.nil? 
      return false, "Remote blast failed"
    else
      if report.iterations.nil? || report.iterations.empty? || report.hits.nil? || report.hits.empty? || report.hits.first.nil?
        return false, "No reciprocal hit by remote blast: no hit in remote database"
      else
        report.hits.each do |hit| # had to add this since a hit may not have a cds marked in the region of the hit, so try hits lower down the e-value ranking
          # puts hit.definition
          if significant_hit?(:hit => hit, :query_sequence => options[:query_sequence], :percent_identity_cutoff => options[:percent_identity_cutoff], :minimum_hit_length => options[:minimum_hit_length], :display_metrics => options[:display_metrics])
            accession = hit.accession
            puts accession
            if options[:embl_objects] && options[:embl_objects].has_key?(accession)
              embl_or_genbank_object = options[:embl_objects][accession]
            else
              embl_or_genbank_object = retrieve_embl_object(accession)
              if embl_or_genbank_object.nil? || embl_or_genbank_object.definition.nil? || embl_or_genbank_object.definition == ""
                embl_or_genbank_object = retrieve_refseq_object(accession) # try retrieving refseq object
                if embl_or_genbank_object.nil? || embl_or_genbank_object.definition.nil? || embl_or_genbank_object.definition == ""
                  next
                end
              end
            end
            if options[:blast_program] == "blastn"
              hit_sequence, hit_qualifiers = retrieve_cds_hit_details_from_rich_sequence_object(hit,embl_or_genbank_object)
            elsif options[:blast_program] == "blastp"
              biosequence = embl_or_genbank_object.to_biosequence
              cds = biosequence.features.select{|feature| feature.feature == "CDS"}.first
              hit_qualifiers = cds.qualifiers
              biosequence = biosequence.auto
              if biosequence.class == Bio::Sequence::NA
                hit_sequence = biosequence.translate(1,11).seq
              elsif biosequence.class == Bio::Sequence::AA
                hit_sequence = biosequence.seq
              end
            end
            next if hit_sequence.nil? # if cds not marked in region of hit try next
            qualifiers_array = hit_qualifiers.map{|qualifier| [qualifier.qualifier, qualifier.value]}
            qualifiers_array.unshift(["organism", "#{embl_or_genbank_object.definition} (#{embl_or_genbank_object.accession})"])
            # do reciprocal blast
            reciprocal_hit_found, number_of_reciprocal_hits_found, message = perform_reciprocal_blast(:query_sequence_for_reciprocal_blast => hit_sequence, :original_query_sequence_id => options[:query_sequence].entry_id, :path_to_blast_executable => options[:path_to_blast_executable], :blast_program => options[:blast_program], :query_sequence_blast_database => options[:query_sequence_blast_database], :blast_options => '-e 1e-20 -F F -b 10 -v 10') 
            if reciprocal_hit_found
              if number_of_reciprocal_hits_found == 1 # just one reciprocal hit to original
                return true, "#{embl_or_genbank_object.definition} (#{embl_or_genbank_object.accession})", qualifiers_array, embl_or_genbank_object
              else # multiple significant reciprocal hits
                return false, "#{embl_or_genbank_object.definition} (#{embl_or_genbank_object.accession})", qualifiers_array, embl_or_genbank_object
              end
            else
              return false, message
            end
          else
            return false, "No reciprocal hit by remote blast: the hit in the remote database was not significant" # no sig hit so not a recipiprocal hit
          end
        end
        return false, "No reciprocal hit by remote blast: hits in remote database but either no hits could be retrieved from EMBL or none contained a CDS in the region of the hit"
      end
    end
  end
  # method to extract the organism qualifiers that were added to the local database entries or the remote hit qualifiers
  # @param [Array] qualifiers An array of two element arrays where the first element in the array is the qualifier type
  # @return [String] Organism description
  def extract_organism_qualifier(qualifiers)
    index_to_delete = nil
    organism = nil
    qualifiers.each_with_index do |qualifier, index|
      qualifier_type = qualifier[0]
      qualifier_value = qualifier[1]
      if qualifier_type == "organism"
        index_to_delete = index
        organism = qualifier_value
        break
      end
    end
    qualifiers.delete_at(index_to_delete)
    return organism
  end
  
  
  # A method to perform the reciprocal blast and see if the reciprocal hits contain the original query sequence
  # @param [Hash] options The options for the method
  # @option options [Bio::Sequence] :query_sequence_for_reciprocal_blast The sequence of the hit which will be used as the query sequence in the reciprocal blast and hopfully find 
  # the ORIGINAL query sequence
  # @option options [String] :path_to_blast_executable The path to the blast executable
  # @option options [String] :blast_program The blast program e.g blastn, blastp
  # @option options [String] :blast_database The path and name of the ORGINAL query blast database
  # @option options [String] :original_query_sequence_id The id (usually accession) of the original query sequence
  # @return [Boolean, Integer, String] If a reciprocal hit was found including the original query sequence, Number of significant reciprocal hits(If more than one significant hit was found - possible gene family?), Message
  def perform_reciprocal_blast(options)
    default_options = {
      :percent_identity_cutoff => 90,
      :minimum_hit_length => 85,
      :blast_options => "-e 1e-20 -F F -b 10 -v 10"
    }
    options.reverse_merge!(default_options)
    reciprocal_report = blast_sequence_locally_without_bioruby(:sequence => options[:query_sequence_for_reciprocal_blast], :path_to_blast_executable => options[:path_to_blast_executable], :blast_program => options[:blast_program], :blast_database => options[:query_sequence_blast_database], :blast_options => options[:blast_options])
     
    if reciprocal_report.nil? ||  reciprocal_report.hits.size == 0
      return false, nil,  "No reciprocal hit by local blast: no hit in query database"
    else
      reciprocal_report_hits = Hash.new
      reciprocal_report.hits.each do |reciprocal_report_hit|
        reciprocal_report_hits[reciprocal_report_hit.accession] = significant_hit?(:hit => reciprocal_report_hit, :query_sequence => options[:query_sequence_for_reciprocal_blast], :percent_identity_cutoff => options[:percent_identity_cutoff], :minimum_hit_length => options[:minimum_hit_length])
      end
      if reciprocal_report_hits[options[:original_query_sequence_id]] # check to see if any of the reciprocal hits match the original
        return true, reciprocal_report_hits.values.grep(true).size, "Matched original sequence"
      else
        return false, nil, "None of the significant reciprocal hits matched the original query sequence"
      end
    end
  end
end

