# -*- coding: utf-8 -*-
require 'parslet'

def make_anchor(s)
  String(s).gsub(' ', '_').gsub(/[^[:alnum:]_]/, '')
end

def escape_string(s)
  # at this point, nothing needs escaped
  # but this provides a hook if needed
  String(s)
end

class RSTTransformer < Parslet::Transform

  # for now, leave the values alone
  rule(:number => simple(:x)) { x }
  rule(:string => simple(:x)) { "'%s'" % [x] }
  rule(:string => sequence(:x)) { "''" }
  rule(:true => simple(:x)) { x }
  rule(:false => simple(:x)) { x }

  # text gets escaped
  rule(:text => simple(:x)) { escape_string(x) }

  # make the anchor
  rule(:level => simple(:l), :anchor => simple(:a)) {
    level = String(l).length
    case level
    when 1
      adornment = '='
    when 2
      adornment = '_'
    when 3
      adornment = '~'
    end
    "\n%s\n\%s" % [a, adornment * a.length]
  }

  # options go back to looking like they did in the fortran file
  # except that trailing comments are dropped
  rule(:option => simple(:o), :value => simple(:v)) { '%s = %s' % [o,v] }
  rule(:option => simple(:o), :value => simple(:v), :trailingcomment => simple(:t)) { '%s = %s' % [o,v] }

  # default values will be syntax highlighted
  rule(:default => sequence(:x)) {
    [".. code-block:: fortran\n", *x].join("\n    ")
  }

  # join comments together on adjacent lines
  rule(:comment => sequence(:x)) {
    c = x.join("\n")
    if c.start_with?("    ")
      "::\n\n" + c
    else
      c
    end
  }

  # put a blank line between sections
  rule(:section => sequence(:x)) { x.join("\n\n") }

  # we completely drop the prelude, as we expect the information that
  # makes sense at the top of the defaults file to be sufficiently
  # different than the information that is at the top of the webpage
  rule(:prelude => sequence(:p), :sections => sequence(:s)) { s.join("\n") }

  # once we make to to the root, just return everything
  rule(:controls => simple(:x)) { x }

end
