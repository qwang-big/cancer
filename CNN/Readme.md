## resize image
```
for i in {1..20};do d=$((100+$i*10));convert ../img.jpg -resize ${d}% $d.jpg;done
for i in {1..5};do d=$((100-$i*10));convert ../img.jpg -resize ${d}% $d.jpg;done
```
## bash
```
sshfs -o reconnect,ServerAliveInterval=15,ServerAliveCountMax=3 192.168.61.8:/prj/ST_PRECISION/OCG/ mnt
```
