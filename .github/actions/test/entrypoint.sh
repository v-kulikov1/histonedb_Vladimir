mysqldatadir=/var/lib/mysql/
#/usr/bin/mysqld_safe
echo "Running mysql in the background"
mysqld_safe --datadir=$mysqldatadir --port=13307 &


cd /var/www/histonedb
echo 'Starting apache2'
apachectl start
cd /var/www
bash db_gen.sh -mysql_db_reinit -histdb_reinit
