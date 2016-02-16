#!/bin/bash

# Copy configuration files
if [ ! -d $HOME/.interimR ] ; then mkdir $HOME/.interimR ; fi
if [ ! -f $HOME/.interimR/user.opts ]       ; then cp ./interimR/cfg/user.opts $HOME/.interimR/user.opts ; fi
if [ ! -f $HOME/.interimR/dfs.datafiles ]   ; then cp ./interimR/cfg/dfs.datafiles $HOME/.interimR/dfs.datafiles ; fi
if [ ! -f $HOME/.interimR/model.variables ] ; then cp ./interimR/cfg/model.variables $HOME/.interimR/model.variables ; fi

# find out original user name and primary group
myusername=$( echo $HOME | sed -e "s;/; ;g" | awk '{ print $NF }' )
mygroup=$( id -g -n $myusername)

# fix permission
chown -R $myusername:$mygroup $HOME/.interimR
chmod -R 755 $HOME/.interimR
