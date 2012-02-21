#!/usr/bin/ruby
## ccs, 06.02.2012
## 
## create the 'setup' files needed for analysis in each folder
require 'fileutils'

@mwNames = ["mw1-round", "mw2-round", "mw3-round", "mw6-round"]
@compareNames = ["mw1", "mw2", "mw3", "mw6"]

## creates a simple setup script that should read like this
## "## variables for the #{name} analysis"
## "##"
## "mwString = #{name}"
## "mwString.Up = toupper(mwString)"
## "## load the master combAna file to do the analysis"
## "source("../combAna.R")"
def createSetupEmu(folder, name)
  script = File.new(File.join(folder, "setup.R"), 'w')
  printComment(script, "variables for the #{name} analysis")
  printComment(script, "")
  printMwString(script,name)
  printComment(script, "load the master combAna file and do the analysis")
  printLine(script, "source(\"../fileNames.R\")\n")
  printLine(script, "source(\"../combAna.R\")\n")
  script.close()
end

## create the setup-plot-implaus script
def createSetupPlotImplaus(folder, name, compareNames)
  script = File.new(File.join(folder, "setup-plot-implaus.R"), 'w')
  printComment(script, "variables for #{name} analysis")
  printComment(script, "compare implaus against other observations")
  printMwString(script, name)
  printCompareString(script, name, compareNames)
  printComment(script, "apply to everything")
  printLine(script, "for(index  in 1:ncompare){\n")
  printCompareNameIndex(script)
  printLine(script, "cat(\"# mwCompare: \" , mwStringCompare, \"\\n\")\n")
  printLine(script, "source(\"../fileNames.R\")\n")
  printLine(script, "source(\"../plotImplaus.R\")\n")
  printLine(script, "}\n\n")
  script.close()
end

## create the setup-grid-implaus script
def createSetupGridImplaus(folder, name, compareNames)
  script = File.new(File.join(folder, "setup-grid-implaus.R"), 'w')
  printComment(script, "variables for #{name} analysis")
  printComment(script, "grid the implaus against other observations")
  printMwString(script, name)
  printCompareString(script, name, compareNames)
  printComment(script, "apply to everything")
  printLine(script, "for(index in 1:ncompare){\n")
  printCompareNameIndex(script)
  printLine(script, "cat(\"# mwCompare: \" , mwStringCompare, \"\\n\")\n")
  printLine(script, "source(\"../fileNames.R\")\n")
  printLine(script, "source(\"../gridCompareImplaus.R\")\n")
  printLine(script, "}\n\n")
  script.close()
end

def printCompareNameIndex(file)
  printLine(file, "mwStringCompare = compareStrings[index]\n")
  printLine(file, "mwStringCompare.Up = toupper(mwStringCompare)\n")
  printLine(file, "\n")
end

def printCompareString(file, name, compareNames)
  printLine(file, "compareStrings <- c(")
  count = 0
  compareNames.each do |name|
    if(count < compareNames.length-1)
      printLine(file, "\"#{name}\",")
    else 
      printLine(file, "\"#{name}\"")
    end
    count = count + 1
  end
  printLine(file, ")\n")
  printLine(file, "ncompare <- length(compareStrings)\n")
end

def printMwString(file, name)
  printLine(file, "mwString = \"#{name}\"\n")
  printLine(file, "mwString.Up = toupper(mwString)\n")
end


def printLine(file, string)
  file.print(string)
end

def printComment(file, string)
  file.print("## #{string}\n")
end

@mwNames.each do |name|
  shortName = name.split("-").first
  p(shortName)
  createSetupEmu(name, shortName )
  createSetupPlotImplaus(name, shortName, @compareNames)
  createSetupGridImplaus(name, shortName, @compareNames)
end
