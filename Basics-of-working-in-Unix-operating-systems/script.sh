#! /bin/bash

if [[ $# -ne 2 || $# -ne 4 ]]; then
    echo "There's a wrong input" >&2
    exit 2
fi
file="/etc/group"
if [[ "$1" == "-f" ]]; then
    shift
    file="$1"
    shift
fi
group=$1

str="$(grep -w "$group" "$file")"

if [[ "$str" == "" ]]; then
   echo "This group does not exist" >&2
   exit 1
fi
if [[ -z $group ]]; then
    echo "There is no matched group" >&2
    exit 2
fi

if ! [[ -f "$file" ]]; then
    echo "File has not been found" >&2
    exit 2
fi
users=`grep "$group": "$file"`
users="${users#*:*:*:}"
current="${users%,*}"
users="${users#*,}"
echo  "$current"
while ! [[ "$users" == "$current" ]]
do
        current="${users%,*}"
        users="${users#*,}"
        echo "$current"
done
exit 0
