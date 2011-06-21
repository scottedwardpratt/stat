# ccs, use this to merge in files from your own local repos
# 1) fetch the repo: git fetch /path/to/repo master:newbranch
# 2) move the stuff into its own subdir (using this script) and for above you'd use
# branch_name=newbranch
# 3) merge the fetched stuff onto the main branch
# git merge newbranch
# git branch -d newbranch # delete the dangling branch
# 
# subdir_name = where you want the your external git tree to go  
# branch_name = the branch you fetched your code into
SUBDIR_NAME=rhicEmu
BRANCH_NAME=madaiAna
git filter-branch -f --index-filter \
	    'git ls-files -s | \
	            sed "s-\t-&'"$SUBDIR_NAME"'/-" | \
	                    GIT_INDEX_FILE=$GIT_INDEX_FILE.new git update-index --index-info && \
	                            mv $GIT_INDEX_FILE.new $GIT_INDEX_FILE
    ' "$BRANCH_NAME"
    
