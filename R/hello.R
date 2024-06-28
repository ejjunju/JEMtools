# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}


#Run devtools::document() to create documentation
#There are three main ways to run roxygen:
# roxygen2::roxygenise()
# devtools::document().
# Ctrl + Shift + D, if you’re using RStudio.


#git
#https://resources.github.com/github-and-rstudio/
#create a respository on git: https://github.com/ejjunju/JEMtools
#copy the url
#configure project as a git project and provide link to git respository

# git status # Allows us to see the status of the files on our branch at any given time. Your file is listed under Untracked files
  # #fatal: detected dubious ownership in repository at 'D:/WORK/Tools/R/Packages/JEMtools'
  #'D:/WORK/Tools/R/Packages/JEMtools' is on a file system that does not record ownership
  #To add an exception for this directory, call:
  #   git config --global --add safe.directory D:/WORK/Tools/R/Packages/JEMtools
# git add . #add files Add your file to the staging area so it’s prepared to become part of the next commit:
# git commit -m "Testing commit on 2024.06.24" # local commit
# git log --oneline #See the history of commits:
# git push #push to git

#to get url to repository
#in aterminal (cmd)
#git config --get remote.origin.url

#sometimes Can't push to GitHub because of large file which I already deleted
# i had to create a new project linked to github url
# edit r.ptoj file to make it a package by adding some lines
# BuildType: Package
# PackageUseDevtools: Yes
# PackageInstallArgs: --no-multiarch --with-keep.source
# PackageRoxygenize: rd,collate,namespace,vignette
