#!/usr/bin/ruby
require 'fileutils'
@mwNames = ["mw1-round", "mw2-round", "mw3-round", "mw6-round"]

@compareNames = ["mw1", "mw2", "mw3", "mw6"]

@subFolders = ["design", "output", "images"]

@rawPath = "#{Dir.pwd}/raw-data/"
@comparePath = "#{Dir.pwd}/raw-data"

@runNameSuffix = "_3par.dat"
@desNameSuffix = "_3par_sorted.dat"

def setupFolder(name)
  top = Dir.pwd
  rawNameSuffix = (name.split("-")[0]).upcase

  p(rawNameSuffix) 

  Dir.chdir(name)
  @compareNames.each {|folder| FileUtils.mkdir(folder)}

  @compareNames.each do |folder| 
    Dir.chdir(folder)
    FileUtils.mkdir(@subFolders) 
    FileUtils.cp(File.join(@rawPath, "design_#{rawNameSuffix}#{@desNameSuffix}"), "./design")
    FileUtils.cp(File.join(@rawPath, "lum_fun_outps_#{rawNameSuffix}#{@runNameSuffix}"), "./output")
    FileUtils.cp(File.join(@rawPath, "metallicity_MV_outputs_#{rawNameSuffix}#{@runNameSuffix}"), "./output")
    folderNameSuffix = folder.upcase
    p(folderNameSuffix)
    FileUtils.cp(File.join(@comparePath, "lum_fun_observed_#{folderNameSuffix}#{@runNameSuffix}"), ".")
    FileUtils.cp(File.join(@comparePath, "metallicity_MV_observed_#{folderNameSuffix}#{@runNameSuffix}"), ".")
    Dir.chdir("..")
  end
  ## and switch back
  Dir.chdir(top)
end

@mwNames.each {|name| FileUtils.mkdir(name) unless File.exists?(name)}

@mwNames.each { |name| setupFolder(name)}
  
  
