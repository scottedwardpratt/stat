#!/usr/bin/ruby
## ccs, cec24@phy.duke.edu
## run @runscript at on each topLevel folder
## if @depth is set to true it is run against each of the inner folders
require 'fileutils'

topLevel = ["mw2-round", "mw3-round", "mw6-round"]
inner = ["mw1", "mw2", "mw3", "mw6"]

@runscript = "../plotRoundVsDesign.R"
@depth = false


def makeImages(mwTop, mwCompare)
  runstring = "R --slave --no-save < #{@runscript}"
  ENV['mwAct'] = mwTop
  ENV['mwActUp'] = mwTop.upcase
  ENV['mwComp'] = mwCompare
  ENV['mwCompUp'] = mwCompare.upcase
  system(runstring)
end

topLevel.each do |topFolder|
  top=Dir.pwd
  Dir.chdir(topFolder)
  topString = topFolder.split("-")[0]

  if(@depth == true)
    inner.each do |folder|
      Dir.chdir(folder)
      FileUtils.cp(File.join(top, @runscript), ".")
      makeImages(topString, folder)
      Dir.chdir("..")
    end
  else 
    makeImages(topString, "")
  end
  Dir.chdir(top)
end


