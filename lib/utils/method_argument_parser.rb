# a class to help parse method arguments and check a) if required arguments are present and b) if they are of the correct type
module MethodArgumentParser
  class Parser

    def initialize(&block)
      @options = []
      instance_eval(&block)
    end

    def option(name, option_opts={})
      option = Option.new(name)
      ## fill in :type
      option_opts[:type] = # normalize
      case option_opts[:type]
      when :boolean, :bool, :flag
        '(True|False)Class'
      when :integer,:int
        '(Fixnum|Integer)'
      when :double,:float
        'Float'
      when :string
        'String'
      when :array
        'Array'
      when :hash
        'Hash'
      when Class
        option_opts[:type].name
      when nil
        nil
      else
        option_opts[:type]
      end


      type_from_default =
      case option_opts[:default]
      when Integer
        'Integer'
      when Numeric
        'Float'
      when TrueClass, FalseClass
        '(True|False)Class'
      when String
        'String'
      when nil
        nil
      else
        option_opts[:default]
      end

      raise ArgumentError, ":type specification and default type don't match (default type is #{type_from_default})" if option_opts[:type] && type_from_default && option_opts[:type] != type_from_default

      option.type = option_opts[:type] || type_from_default
      option.default = option_opts[:default] || nil
      if option.default.nil?
        option.required = option_opts[:required] || false
      else
        option.required = true
      end
      @options << option
    end

    def parse(options = {})
      final_options = Hash.new
      @options.each do |option|
        final_options[option.name] = option.default unless option.default.nil?
        unless options[option.name].nil?
          if option.type
            raise ArgumentError, "the option #{option.name} is of the wrong type it should be of type #{option.type}" unless options[option.name].class.name.match(/#{option.type}/)
          end
          final_options[option.name] = options[option.name]
        end
        raise ArgumentError, "the required argument #{option.name} has not been specified" if option.required && !final_options.keys.include?(option.name)
      end
      final_options.merge!(options)
    end

    def self.check_options(options_to_check, &block)
      parser = MethodArgumentParser::Parser.new(&block)
      parser.parse(options_to_check)
    end

    class Option
      attr_accessor :name, :type, :required, :default
      def initialize(name)
        @name = name
      end
    end
  end
end