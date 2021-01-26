
class Command(BaseCommand):
    help = 'Creating histone types and variants for the seed sequences'

    # Logging info
    logging.basicConfig(filename='log/buildvarianttypes.log',
                        format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    log = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument(
            "-f",
            "--force",
            default=False,
            action="store_true",
            help="Force the regeneration of HMM from seeds, HUMMER search in db_file, Test models and loading of results to database")

        parser.add_argument(
            "--profile",
            default=False,
            action="store_true",
            help="Profile the command")

    def _handle(self, *args, **options):
        self.log.info('=======================================================')
        self.log.info('===          buildvarianttypes START           ===')
        self.log.info('=======================================================')
        self.start_time=datetime.now()

        if options["force"]:
            #Clean the DB, removing all types and variants
            Variant.objects.all().delete()
            Histone.objects.all().delete()

        self.log.info('=======================================================')
        self.log.info('===   buildvariants_parallel SUCCESSFULLY finished  ===')
        self.log.info('=======================================================')

    def handle(self, *args, **options):
        if options['profile']:
            profiler = Profile()
            profiler.runcall(self._handle, *args, **options)
            profiler.print_stats()
        else:
            self._handle(*args, **options)
