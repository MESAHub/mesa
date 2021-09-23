# -*- coding: utf-8 -*-

require_relative 'defaultsParser.rb'
require_relative 'rstTransformer.rb'

DEFAULTS_FILES=["star_job.defaults",
                "controls.defaults",
                "pgstar.defaults",
                "binary_controls.defaults",
                ]

OUTPUT_PATH = "./"

def convert_file(filename)
  contents = File.open(filename).read
  begin
    tree = DefaultsParser.new.parse(contents)
  rescue Parslet::ParseFailed => error
    puts error.parse_failure_cause.ascii_tree
  end
  RSTTransformer.new.apply(tree)
end

def write_file(filename, contents)
  f = File.open(filename, "w")
  f.write(contents)
  f.close
end

def generate_filename(filename)
  File.basename(filename).chomp(".defaults") + ".rst"
end

# convert docs
DEFAULTS_FILES.each do |defaults_filename|

  puts "Converting %s ..." % [defaults_filename]

  rst = convert_file(defaults_filename)
  rst_filename = generate_filename(defaults_filename)

  write_file(File.join(OUTPUT_PATH, rst_filename), rst)
  
end
