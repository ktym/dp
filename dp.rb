#!/usr/bin/env ruby
#
# Usage:
#  % ruby dp.rb
#  % ruby dp.rb 'ACCAGT' 'ACAGC'
#

require 'pp'

# Generic interface
class DP

  attr_accessor :gap_penalty, :ext_penalty, :match, :mis_match

  def initialize(target, query)
    @gap_penalty = 2
    @ext_penalty = 1
    @match = 1
    @mis_match = -1
    @target = Sequence.new(target)
    @query = Sequence.new(query)
  end

  def s(x, y)
    return x == y ? @match : @mis_match
  end

  def init_matrix
    raise NotImplementedError.new("#{self.class}##{__method__}: Create @matrix = Matrix.new(tlen, qlen)")
  end

  def calc_matrix
    raise NotImplementedError.new("#{self.class}##{__method__}: Update @matrix with .get() and .set()")
  end

  def trace_back
    raise NotImplementedError.new("#{self.class}##{__method__}: Generate @target_align, @query_align")
  end

  class Sequence
    def initialize(seq)
      @seq = seq
    end

    def to_s
      @seq
    end

    def [](one_base_position)
      zero_base_position = one_base_position - 1
      @seq[zero_base_position]
    end

    def length
      @seq.length
    end
  end

  class Matrix
    def initialize(x, y, default = 0)
      @matrix = Array.new(x+1) { Array.new(y+1, default) }
    end

    def get(i, j)
      @matrix[j][i]
    end

    def set(i, j, v)
      @matrix[j][i] = v
    end
  end

  def dp
    puts "Target: #{@target}"
    puts "Query:  #{@query}"
    init_matrix
    pp @matrix if $DEBUG
    calc_matrix
    pp @matrix
    trace_back
    puts "Target: #{@target_align}"
    puts "Query:  #{@query_align}"
  end

end

# Global alignment
class NeedlemanWunsch < DP

  def init_matrix
    @matrix = Matrix.new(@query.length, @target.length, 0)
    1.upto(@target.length) do |i|
      value = - (i * @gap_penalty)
      @matrix.set(i, 0, value)
    end
    1.upto(@query.length) do |j|
      value = - (j * @gap_penalty)
      @matrix.set(0, j, value)
    end
  end

  def calc_matrix
    1.upto(@target.length) do |i|
      1.upto(@query.length) do |j|
        value = [
          @matrix.get(i-1, j-1) + s(@target[i], @query[j]),
          @matrix.get(i,   j-1) - @gap_penalty,
          @matrix.get(i-1, j  ) - @gap_penalty,
        ].max
        @matrix.set(i, j, value)
      end
    end
  end

  def trace_back
    i = @target.length
    j = @query.length
    target_trace = ""
    query_trace = ""

    score = @matrix.get(i, j)
    while i > 0 or j > 0
      case score
      when @matrix.get(i-1, j) - @gap_penalty
        score = @matrix.get(i-1, j)
        target_trace += @target[i]
        query_trace += '-'
        i -= 1
      when @matrix.get(i, j-1) - @gap_penalty
        score = @matrix.get(i, j-1)
        target_trace += '-'
        query_trace += @query[j]
        j -= 1
      else
        score = @matrix.get(i-1, j-1)
        target_trace += @target[i]
        query_trace += @query[j]
        i -= 1
        j -= 1
      end
    end

    @target_align = target_trace.reverse
    @query_align = query_trace.reverse
  end

end

# Local alignment
class SmithWaterman < DP

  def init_matrix
    @matrix = Matrix.new(@query.length, @target.length, 0)
  end

  def calc_matrix
    @best = {:score => 0, :i => 0, :j => 0}
    1.upto(@target.length) do |i|
      1.upto(@query.length) do |j|
        value = [
          @matrix.get(i-1, j-1) + s(@target[i], @query[j]),
          @matrix.get(i,   j-1) - @gap_penalty,
          @matrix.get(i-1, j  ) - @gap_penalty,
          0,
        ].max
        @matrix.set(i, j, value)
        if @matrix.get(i, j) > @best[:score]
          @best[:score] = @matrix.get(i, j)
          @best[:i] = i
          @best[:j] = j
        end
      end
    end
  end

  def trace_back
    i = @best[:i]
    j = @best[:j]
    target_trace = ""
    query_trace = ""

    score = @matrix.get(i, j)
    while i > 0 or j > 0
      case score
      when @matrix.get(i-1, j) - @gap_penalty
        score = @matrix.get(i-1, j)
        target_trace += @target[i]
        query_trace += '-'
        i -= 1
      when @matrix.get(i-1, j) - @gap_penalty
        score = @matrix.get(i, j-1)
        target_trace += '-'
        query_trace += @query[j]
        j -= 1
      else
        score = @matrix.get(i-1, j-1)
        target_trace += @target[i]
        query_trace += @query[j]
        break if score == 0
        i -= 1
        j -= 1
      end
    end

    @target_align = target_trace.reverse
    @query_align = query_trace.reverse
  end

end

# Global alignment with affine gap penalty
class NeedlemanWunschGotoh < DP

  Inf = Float::INFINITY

  def init_matrix
    @matrix = {
      :m => Matrix.new(@query.length, @target.length),
      :x => Matrix.new(@query.length, @target.length),
      :y => Matrix.new(@query.length, @target.length),
    }
    @matrix[:m].set(0, 0, 0)
    @matrix[:x].set(0, 0, -Inf)
    @matrix[:y].set(0, 0, -Inf)
    1.upto(@target.length) do |i|
      @matrix[:m].set(i, 0, -Inf)
      @matrix[:y].set(i, 0, -Inf)
      value = - @gap_penalty - @ext_penalty * (i-1)
      @matrix[:x].set(i, 0, value)
    end
    1.upto(@query.length) do |j|
      @matrix[:m].set(0, j, -Inf)
      @matrix[:x].set(0, j, -Inf)
      value = - @gap_penalty - @ext_penalty * (j-1)
      @matrix[:y].set(0, j, value)
    end
  end

  def calc_matrix
    1.upto(@target.length) do |i|
      1.upto(@query.length) do |j|
        s_value = s(@target[i], @query[j])
        value_m = [
          @matrix[:m].get(i-1, j-1) + s_value,
          @matrix[:x].get(i-1, j-1) + s_value,
          @matrix[:y].get(i-1, j-1) + s_value,
        ].max
        @matrix[:m].set(i, j, value_m)

        value_x = [
          @matrix[:m].get(i-1, j) - @gap_penalty,
          @matrix[:x].get(i-1, j) - @ext_penalty,
        ].max
        @matrix[:x].set(i, j, value_x)

        value_y = [
          @matrix[:m].get(i, j-1) - @gap_penalty,
          @matrix[:y].get(i, j-1) - @ext_penalty,
        ].max
        @matrix[:y].set(i, j, value_y)
      end
    end
  end

  def trace_back
    i = @target.length
    j = @query.length

    score = {
      :m => @matrix[:m].get(i, j),
      :x => @matrix[:x].get(i, j),
      :y => @matrix[:y].get(i, j),
    }
    cur = score.max{|a,b| a[1] <=> b[1]}.first
    target_trace = ""
    query_trace = ""

    while i > 0 or j > 0
      score = @matrix[cur].get(i, j)
      s_value = s(@target[i], @query[j])
      case cur
      when :m
        case score
        when @matrix[:m].get(i-1, j-1) + s_value
          score = @matrix[:m].get(i-1, j-1)
          cur = :m
        when @matrix[:x].get(i-1, j-1) + s_value
          score = @matrix[:x].get(i-1, j-1)
          cur = :x
        when @matrix[:y].get(i-1, j-1) + s_value
          score = @matrix[:y].get(i-1, j-1)
          cur = :y
        end
        target_trace += @target[i]
        query_trace += @query[j]
        i -= 1
        j -= 1
      when :x
        case score
        when @matrix[:m].get(i-1, j) - @gap_penalty
          score = @matrix[:m].get(i-1, j)
          cur = :m
        when @matrix[:x].get(i-1, j) - @ext_penalty
          score = @matrix[:x].get(i-1, j)
          cur = :x
        end
        target_trace += @target[i]
        query_trace += '-'
        i -= 1
      when :y
        case score
        when @matrix[:m].get(i, j-1) - @gap_penalty
          score = @matrix[:m].get(i, j-1)
          cur = :m
        when @matrix[:y].get(i, j-1) - @ext_penalty
          score = @matrix[:y].get(i, j-1)
          cur = :y
        end
        target_trace += '-'
        query_trace += @query[j]
        j -= 1
      end
    end

    @target_align = target_trace.reverse
    @query_align = query_trace.reverse
  end

end

if __FILE__ == $0

  target = ARGV.shift || 'ACCAGT'
  query = ARGV.shift || 'ACAGC'

  puts "# Needleman-Wunsch"
  nw = NeedlemanWunsch.new(target, query)
  nw.dp

  puts "# Smith-Waterman"
  sw = SmithWaterman.new(target, query)
  sw.dp

  puts "# Needleman-Wunsch-Gotoh"
  nwg = NeedlemanWunschGotoh.new(target, query)
  nwg.dp

end


