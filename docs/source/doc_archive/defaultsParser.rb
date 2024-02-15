# -*- coding: utf-8 -*-
require 'parslet'

class DefaultsParser < Parslet::Parser

  # the format being parsed is described in 
  # $MESA_DIR/star/defaults/FORMAT

  # simple rules that we'll reuse throughout

  rule(:equals)     { space? >> str('=') >> space? }

  rule(:bang)       { str('!') }
  rule(:hash)       { str('#') }
  rule(:quote)      { str('"') | str("'") }

  rule(:space1)     { match("[ \t]") }
  rule(:space)      { space1.repeat(1) }
  rule(:space?)     { space.maybe }

  rule(:newline)    { str("\n") }
  rule(:empty)      { space? >> newline }
  rule(:empty?)     { empty.maybe }


  # rules for parsing a default
  # these look like option = value
  # option is 
  # a value is a string, number, or boolean

  rule(:option)     { match('[A-Za-z0-9_(:),]').repeat(1).as(:option) }

  rule(:digit)      { match('[0-9]') }

  rule(:number) {
    (str('-').maybe >>
      (str('0') | (match('[1-9]') >> digit.repeat)) >>
      (str('.') >> digit.repeat(1)).maybe >>
      (match('[eEdD]') >> (str('+') | str('-')).maybe >> digit.repeat(1)).maybe
     ).as(:number)
  }

  rule(:string) {
    quote >> (str('\\') >> any | quote.absent? >> any).repeat.as(:string) >> quote
  }


  rule(:value) {
    string | number |
    str('.true.').as(:true) | str('.false.').as(:false)
  }

  rule(:default) { space? >> option >> equals >> value.as(:value) >> (comment.as(:trailingcomment) | empty) }


  # a divider looks like !---- (or more)
  rule(:divider) { bang >> match('-').repeat(4) }

  # an anchor looks like !# some text
  rule(:anchor) { space? >> bang >> hash.repeat(1,3).as(:level) >> space1 >> text.as(:anchor) }

  # a regular comment looks like ! some text
  rule(:text) { match('[^\r\n]').repeat(1) }
  rule(:comment) { space? >> bang >> space1 >> text.as(:text) >> empty }

  # an empty comment looks like
  rule(:emptycomment) { space? >> bang >> space1.maybe >> newline }

  # the whole file comes in blocks

  # a block can be:
  #   a single anchor
  #   multiple consecutive comment lines
  #   a single defaults line 

  rule(:block) { empty.repeat >> ( anchor |
                                   comment.repeat(1).as(:comment) |
                                   default.repeat(1).as(:default) ) >> (emptycomment | empty).repeat }

  # a section is composed of multiple blocks
  rule(:section) { divider >> block.repeat(1).as(:section) >> empty? }

  # there is a special section at the start called the prelude which
  # is all comments or empty lines
  rule(:preludeline) { comment | empty }
  rule(:prelude) { preludeline.repeat(1) }

  rule(:controls) { prelude.as(:prelude) >> section.repeat(1).as(:sections) }
  root :controls

end
