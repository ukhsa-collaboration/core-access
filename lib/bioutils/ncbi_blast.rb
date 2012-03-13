module Bio
  class Blast
    # method to perform a remote blast on the ncbi site
    class NCBI
      require 'utils/hash_reverse_merge'

      # create a module for catching Net::Http errors
      module Net::HTTPBroken
      end

      # Include the module into the exceptions
      [Timeout::Error, Errno::EINVAL, Errno::ECONNRESET, Errno::ECONNREFUSED, EOFError, Net::HTTPBadResponse, Net::HTTPHeaderSyntaxError, Net::ProtocolError, SocketError].each {|m| m.send(:include, Net::HTTPBroken)}


      def initialize(options)
        @host = "www.ncbi.nlm.nih.gov"
        @path = "/blast/Blast.cgi"
        options[:submit_params][:program]

        default_options = {
          :blast_program => 'blastn',
          :blast_database => 'nr',
          :submit_params => {},
          :retrieve_params => {},
          :retries => 120 # will retry each net::http operation for a default of 120 times (with the sleep of 2 this results in a minute of retries)
        }
        options.replace(default_options.merge(options)) # reverse merge defaults into options
        options[:matrix_name] = 'blosum62' if (options[:blast_program] == 'blastp' || options[:blast_program] == 'blastx' || options[:blast_program] == 'tblastx') && (options[:matrix_name].nil? || options[:matrix_name] == "")
        @retries = options[:retries]
        options[:submit_params][:program] = options[:blast_program]
        options[:submit_params][:database] = options[:blast_database]
        # These parameters are specified in this doc http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.pdf
        default_submit_params = {
          'CMD'                           => 'Put',
          'FORMAT_OBJECT'                 => 'Alignment',
          'COMPOSITION_BASED_STATISTICS'  => 'off', 
          'DATABASE'                      => '',
          'EXPECT'                        => '1e-10', 
          'FILTER'                        => '',
          'PROGRAM'                       =>  '',
          'SERVICE'                       => 'plain',
          'QUERY'                         => '',
          'OTHER_ADVANCED'                => ''
        }
        default_retrieve_params = {
          'CMD'            => 'Get',
          'DESCRIPTIONS'   => 100,
          'FORMAT_TYPE'    => 'XML',
          'ALIGNMENTS'     => 50,
          'ALIGNMENT_VIEW' => 'Tabular'
        }

        # parse options and merge with defaults
        @submit_params = parse_params_hash_and_merge_with_defaults(options[:submit_params], default_submit_params)
        @retrieve_params = parse_params_hash_and_merge_with_defaults(options[:retrieve_params], default_retrieve_params)

      end
      def query(query_sequence)
        case query_sequence
        when Bio::Sequence
          query_sequence = query.output(:fasta)
        when Bio::Sequence::NA, Bio::Sequence::AA, Bio::Sequence::Generic
          query_sequence = query_sequence.to_fasta('query_sequence', 70)
        else
          query_sequence = query_sequence.to_s
        end

        @submit_params['QUERY'] = CGI.escape(query_sequence)

        data = []

        @submit_params.each do |k, v|
          data.push("#{k}=#{v}") if v
        end

        report = nil

        # begin
        success = false
        while !success
          http = Bio::Command.new_http(@host)
          http.open_timeout = 300
          http.read_timeout = 600
          retries = @retries
          begin
            result, = http.post(@path, data.join('&'))
          rescue Net::HTTPBroken => e
            if retries > 0
              retries -=1
              sleep 1 and retry
            else
              puts e.message
            end
          end
          output = result.body
          # find RID
          rid =""
          output.each_line do |line|
            if line =~ /RID = (\w+)\n$/
              rid = $1
              break
            end
          end
          retrieve_param_data = []

          @retrieve_params.each do |k, v|
            retrieve_param_data.push("#{k}=#{v}") if v
          end
          # push on RID value
          retrieve_param_data.push("RID=#{rid}")
          status = "WAITING"
          # retrieve result
          until status !~ /WAITING/
            sleep 5
            begin
              result, = http.post(@path, retrieve_param_data.join('&'))
            rescue Net::HTTPBroken => e
              if retries > 0
                retries -=1
                sleep 1 and retry
              else
                puts e.message
              end
            end
            output = result.body
            status = ""
            output.each_line do |line|
              if line =~ /Status=(.+)\n$/
                status = $1
                break
              elsif line =~ /<\?xml version=/ # account for XML format returns
                status = "XML_READY"
                break
              end
            end
            status = "NO STATUS" if status == ""
          end
          # if READY reponse parse result and return report
          if status == "READY"
            raw_txt = output.split(/\<\/?PRE\>/)[1]
            text_for_blast_report =""
            raw_txt.each_line do |line|
              text_for_blast_report += line unless line =~ /^#/ or line !~ /\t/# remove blank and commented lines
            end
            success = true
          elsif status == "XML_READY"
            text_for_blast_report = output
            success = true
          else
            puts "No result was obtained: status was #{status}"
            return nil
          end
        end
        return Bio::Blast::Report.new(text_for_blast_report)
      end

      private
      def parse_params_hash_and_merge_with_defaults(params, default_params)
        converted_params = Hash.new
        params.each_pair do |key, value|
          if key.is_a?(Symbol)
            params.delete(key)
            key = key.to_s.upcase
            converted_params[key] = value
          end
        end
        #end
        params.merge!(converted_params)
        params.replace(default_params.merge(params)) # reverse_merge! hash
        return params
      end

    end
  end
end
