#!/usr/bin/ruby
## ccs, cec24@phy.duke.edu
## 
## take a list of the variables to analyse and then create a comb.dat file for them. 
## 
## Essentially, when we get the raw data it is stored in the rows of many disparate files, one in each run 
## folder. Where a run-folder corresponds to a location in the design space. All we do here is to collect
## the rows we want from all the run-folders and present them as columns in a final combined output file.
## 
## Each row in the final output gives the specified observables for a given run. The file is sorted in 
## run order, the first also gives the run-id (this can be turned off here)
##
##
## input: 
## - path to process
## - centrality classes
##   cclasses go in steps of 10, 1 is 0-10, 2 is 10-20, 3 is 20-30 etc
## - list of variable names to collect (from stdin), check testSpec.dat for an example
##
## output:
## - combined data file
## - an errors file can be produced as an option
## - the default output file
## the vertical fields in this file will be given in the same order as they are in the specArray

@printRunId = false


## holds the observables we care about for a given run.
## runFolder -> path to the results.txt file we care about, this therefore includes the centclass path
## runIndex -> run id, needed for sorting the runs
## specArray -> the field names of the output we're going to keep
class RunResultsHolder
  attr_reader :results, :index 
  def initialize(runFolder, runIndex, specArray)
    @index = runIndex.to_i
    @dataFile = "results.dat" ## change this if you need to process a different kind of file
    @results = scrapeResults(runFolder, specArray)
  end
  ## reads the results.dat file in folder and stores the lines
  ## matching the values in specArray
  def scrapeResults(folder, specArray)
    tempResults = []
    File.open(File.join(folder, @dataFile)).each do |line|
      # split the line on spaces, the 2nd field should match one of the 
      # fields in specArray if can use it
      fieldName = line.split(" ")[1]
      if( specArray.index(fieldName) != nil ) 
        rtemp = {}
        rtemp[:name] = fieldName
        # we should track what order in the spec this particular result is
        rtemp[:index] = specArray.index(fieldName)
        rtemp[:value] = line.split(" ")[2]
        rtemp[:err] = line.split(" ")[3]
        tempResults.push(rtemp)
      end
    end
    ## now we sort the results by their index
    tempResults.sort { |x,y| x[:index] <=> y[:index]}
  end

  def printResultField(file, fieldKey)
    if(fieldKey == :value || fieldKey == :err)
      if(@printRunId)
        file.print "#{@index} "
      end
    end
    @results.each {|res| file.print "#{res[fieldKey]} " }
    file.print "\n"
  end
      
  def printResults(outfile)
    printResultField(outfile, :value)
  end
  
  def printResultErrors(outfile)
    printResultField(outfile, :err)
  end

  def printHeader(outfile)
    printResultField(outfile, :name)
  end
    
end

## read which the observable specification from stdin
def readDataSpecs(specfile)
  centClasses = ["cent0to5_","cent20to30_"] #These are the centralities used by allb
  specs = []
  if(ARGV[1]=="allb")
    STDIN.each do |line|
      centClasses.each do |centprefix|
        linetemp=centprefix+line
        if(linetemp.match(/^#\w*/) == nil && linetemp.match(/^s/) == nil)
          specs.push(linetemp.chomp!)
        end
      end
    end
  else
    STDIN.each do |line|
      if(line.match(/^#\w*/) == nil && line.match(/^s/) == nil)
        specs.push(line.chomp!)
      end
    end
  end

  if(specs == [])
    stderr.print("# error didn't read any valid field names ")
    exit -1
  end
  specs
end

## print some header info to the outputfile
def outputHeader(outputfile, centFolderString)
  outputfile.print("# combined data from runs in #{@processPath}\n")
  outputfile.print("# centrality from #{centFolderString}\n#\n")
  outputfile.print("# generated by: generate-analysis-files.rb on #{Time.now}\n")
  outputfile.print("#\n")
end

## print some header info to the outputfile
def outputDefaultHeader(outputfile, centFolderString)
  outputfile.print("# default data for runs in #{@processPath}\n")
  outputfile.print("# centrality from #{centFolderString}\n#\n")
  outputfile.print("# generated by: generate-analysis-files.rb on #{Time.now}\n")
  outputfile.print("#\n# first line is the value, second line gives  the errors\n")
  outputfile.print("#\n")
end



## ## ## ## ## ## ## ## ## ## ## 
## RUNTIME CODE STARTS  HERE  ##
## ## ## ## ## ## ## ## ## ## ##

## user feedback
if(ARGV.length < 2) 
  puts "# generate-analysis-files"
  puts "# if no outputfile given, then will use stdout"
  puts "# if errfile arg is given an outputfile with errors instead of values will be given"
  puts "# arguments <process-path> <centrality-class> [<outputfile> <errfile>]"
  exit -1
end

## process cmdline args, error output is optional
@processPath = File.expand_path(ARGV[0])
@centFolderString = ARGV[1] 
if(ARGV.length > 2)
  @finalOutputFile = ARGV[2]
  if(ARGV.length > 3)
    @errorOutputFile = ARGV[3]
  end
else 
  @finalOutputFile = STDOUT
  @errorOutputFile = nil
end

unless File.exists?(@processPath) && File.directory?(@processPath)
  puts "#{@processPath} doesn't exist or isn't a dir"
  exit -1
end

puts "# collating data from: #{@processPath}"
puts "# in cent range: #{@centFolderString}"
puts "# results output in: #{@finalOutputFile}"
if(@errorOutputFile)  
  puts "# errors output in: #{@errorOutputFile}" 
end

puts "# reading specification from STDIN"
@dataSpec = readDataSpecs(STDIN)

top=Dir.pwd 
Dir.chdir(@processPath)
## loop over runs, for each run we create a RunResultsHolder object
allResults = []
Dir.glob("run*").each do |runfolder| 
  ## extract the run index
  ## the first entry is the match string, the second is the first match to our () etc
  index = runfolder.match(/run([0-9]*)/)[1]
  #p([runfolder, index, Dir.pwd])
  rtemp = RunResultsHolder.new(File.join(runfolder, @centFolderString),
   index, @dataSpec)
  allResults.push(rtemp) unless (rtemp.results == []) 
end

if(allResults == [])
  puts "Reading data failed!?!?!?"
end

## sort results on their index, otherwise it'll be impossible to match up against
## the design
allResults = allResults.sort{|x,y| x.index <=> y.index }

## now collect the default information
defaultResults = []
defaultFolder = "./default"
defaultTemp = RunResultsHolder.new(File.join(defaultFolder, @centFolderString), 
                                   1, @dataSpec)
defaultResults.push(defaultTemp) unless (defaultTemp.results == [])



## cd back to where we started, it's confusing for the user if the 
## output files don't appear relative to their path
Dir.chdir(top)

if(@finalOutputFile != STDOUT)
  outputFile = File.open(@finalOutputFile, "w")
else
  outputFile = STDOUT
end

outputHeader(outputFile, @centFolderString)
allResults[0].printHeader(outputFile)
allResults.each {|res| res.printResults(outputFile)}

## print the default information
#defaultOutputFileName = @finalOutputFile.split(".").shift +  "-default.dat"
defaultOutputPath = File.dirname(@finalOutputFile)
defaultOutputFile = File.basename(@finalOutputFile).split(".").shift + "-default.dat"
defaultOutputFileName = File.join(defaultOutputPath, defaultOutputFile)


puts("# default data to: #{defaultOutputFileName}")
defaultFile = File.open(defaultOutputFileName, "w")
outputDefaultHeader(defaultFile, @centFolderString)
defaultResults[0].printHeader(defaultFile)
defaultResults[0].printResults(defaultFile)
defaultResults[0].printResultErrors(defaultFile)  
defaultFile.close()


if(outputFile != STDOUT)
  outputFile.close()
end


## now handle error output if any
if(@errorOutputFile != nil)
  outputFile = File.open(@errorOutputFile, "w")
  outputHeader(outputFile, @centFolderString)
  outputFile.print("# errors for given obs \n" )
  allResults[0].printHeader(outputFile)
  allResults.each {|res| res.printResultErrors(outputFile)}
  outputFile.close()
end
                   
  
