from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from browse.models import HistoneDescription, Histone, Variant, OldStyleVariant, Publication
from djangophylocore.models import Taxonomy

import json, os, logging, configparser
from datetime import datetime
from bibtexparser.bparser import BibTexParser

config = configparser.ConfigParser()
config.read('./histonedb.ini')

class Command(BaseCommand):
    help = 'Creating histone types and variants for the seed sequences'

    # Logging info
    logging.basicConfig(filename=os.path.join(config['LOG']['database_log'], "buildhistoneclasses.log"),
                        format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    log = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument(
            "--profile",
            default=False,
            action="store_true",
            help="Profile the command")

    def _handle(self, *args, **options):
        self.log.info('=======================================================')
        self.log.info('===           buildhistoneclasses START             ===')
        self.log.info('=======================================================')
        self.start_time=datetime.now()

        #Clean the DB, removing all types and variants
        Histone.objects.all().delete()
        Variant.objects.all().delete()
        HistoneDescription.objects.all().delete()
        Publication.objects.all().delete()

        with open(config['DATA']['variants'], encoding='utf-8') as f:
            self.variants_json = json.load(f)

        self.create_publications()

        # Populate our Histone variants table add descriptions
        self.create_histone_classes()

        self.log.info('=======================================================')
        self.log.info('===    buildhistoneclasses SUCCESSFULLY finished    ===')
        self.log.info('=======================================================')

    def handle(self, *args, **options):
        if options['profile']:
            from cProfile import Profile
            profiler = Profile()
            profiler.runcall(self._handle, *args, **options)
            profiler.print_stats()
        else:
            self._handle(*args, **options)

    def create_publications(self):
        with open(config['DATA']['publications'], encoding='utf-8') as f:
            bibdb = BibTexParser().parse_file(f)
        for ent in bibdb.entries:
            try:
                publication = Publication.objects.create(id=ent.get('ID'),
                                                         title=ent.get('title'),
                                                         volume=ent.get('volume'),
                                                         issn=ent.get('issn'),
                                                         url=ent.get('url'),
                                                         doi=ent.get('doi'),
                                                         pages=ent.get('pages'),
                                                         number=ent.get('number'),
                                                         journaltitle=ent.get('journaltitle'),
                                                         author=ent.get('author'),
                                                         date=ent.get('date'),
                                                         cited=False)
                self.log.info(f"Publication {publication.id} was created in database.")
            except:
                self.log.error(f"Publication {ent.get('ID')} was not created in database.")
                raise

    def create_description(self, hclass):
        hclass_desc = self.variants_json['info'][hclass]['description']
        try:  # this is temporary hack for old style description, this try-except should be removed after approving classification.json
            obj = HistoneDescription.objects.create(summary=hclass_desc['summary'],
                                                    taxonomy=hclass_desc['taxonomy'],
                                                    genes=hclass_desc['genes'],
                                                    evolution=hclass_desc['evolution'],
                                                    expression=hclass_desc['expression'],
                                                    knock_out=hclass_desc['knock-out'],
                                                    function=hclass_desc['function'],
                                                    sequence=hclass_desc['sequence'],
                                                    localization=hclass_desc['localization'],
                                                    deposition=hclass_desc['deposition'],
                                                    structure=hclass_desc['structure'],
                                                    interactions=hclass_desc['interactions'],
                                                    disease=hclass_desc['disease'],
                                                    caveats=hclass_desc['caveats'])
        except TypeError:
            obj = HistoneDescription.objects.create(summary=hclass_desc)
        self.log.info(f"Description of histone {hclass} type was created in database with id {obj.id}.")
        return obj

    def add_publications(self, hclass_obj):
        for publication_id in self.variants_json['info'][hclass_obj.id]["publications"]:
            try:
                publication = Publication.objects.get(id=publication_id)
                hclass_obj.publications.add(publication)
            except:
                self.log.error(publication_id)
                raise

    def create_histone_classes(self):
        """Create basic histone types"""
        for htype in self.variants_json['tree'].keys():
            obj = self.create_description(hclass=htype)
            obj = Histone.objects.create(id=htype, taxonomic_span=self.variants_json['info'][htype]['taxonomic_span'],
                                         description=obj)
            self.log.info("Histone {} type was created in database.".format(obj.id))

            # self.add_publications(hclass_obj=obj)

            self.create_histone_variants(variants=self.variants_json['tree'][htype], hist_type=htype)

    def create_histone_variants(self, variants, hist_type, parent=None):
        """Create variants (including generics for each histone type) listed in variants_list.json"""
        if variants=="null": return
        for variant in variants.keys():
            obj = self.create_description(hclass=variant)
            obj = Variant.objects.create(id=variant, hist_type_id=hist_type,
                                         taxonomic_span=self.variants_json['info'][variant]['taxonomic_span'],
                                         doublet=False, description=obj, parent_id=parent)
            self.log.info("Created {} variant model in database".format(obj.id))

            # self.add_publications(hclass_obj=obj)

            for alternate_name in self.variants_json['info'][variant]["alternate_names"]:
                tax_name = alternate_name.get("taxonomy", "eukaryota")
                # tax_eu = Taxonomy.objects.get(name='eucaryotes')
                # print('===================== {}'.format(tax_eu))
                alt_variant = OldStyleVariant(
                    updated_variant=obj,
                    name=alternate_name["name"],
                    gene=alternate_name.get("gene"),
                    splice=alternate_name.get("splice"),
                    taxonomy=Taxonomy.objects.get(name=tax_name)
                )
                try:
                    alt_variant.save()
                except:
                    from django.db import connection
                    cursor = connection.cursor()
                    cursor.execute(
                        "ALTER TABLE browse_feature CONVERT TO CHARACTER SET utf8 COLLATE     utf8_general_ci;")
                    alt_variant.save()

            self.create_histone_variants(variants=variants[variant], hist_type=hist_type, parent=variant)
