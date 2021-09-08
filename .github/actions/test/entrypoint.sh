mysqldatadir=/var/lib/mysql/
#/usr/bin/mysqld_safe
echo "Running mysql in the background"
mysqld_safe --datadir=$mysqldatadir --port=13307 &

rm -rf /var/www/histonedb
ln -s /github/workspace /var/www/histonedb


cd /var/www/histonedb
echo 'Starting apache2'
apachectl start
cd /var/www
bash -e db_gen.sh -mysql_db_reinit -histdb_reinit
