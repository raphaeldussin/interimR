#!/bin/bash

# Copy configuration files
SRCDIR=$( pwd )
if [ ! -d $HOME/.interimR ] ; then mkdir $HOME/.interimR ; else echo "directory $HOME/.interimR already exists" ; fi
if [ ! -f $HOME/.interimR/user.opts ]       ; then cat ./interimR/cfg/user.opts | sed -e "s;<SRCDIR>;$SRCDIR;g" \
                                                       > $HOME/.interimR/user.opts ; else 
   echo "file $HOME/.interimR/user.opts already exists, skipping" ;fi
if [ ! -f $HOME/.interimR/dfs.datafiles ]   ; then cat ./interimR/cfg/dfs.datafiles | sed -e "s;<SRCDIR>;$SRCDIR;g" \
                                                       > $HOME/.interimR/dfs.datafiles ; else

   echo "file $HOME/.interimR/dfs.datafiles already exists, skipping" ;fi
if [ ! -f $HOME/.interimR/model.variables ] ; then cp ./interimR/cfg/model.variables $HOME/.interimR/model.variables ; else
   echo "file $HOME/.interimR/model.variables already exists, skipping" ;fi

# find out original user name and primary group
myusername=$( echo $HOME | sed -e "s;/; ;g" | awk '{ print $NF }' )
mygroup=$( id -g -n $myusername)

# fix permission
chown -R $myusername:$mygroup $HOME/.interimR
chmod -R 755 $HOME/.interimR
