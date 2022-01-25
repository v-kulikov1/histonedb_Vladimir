from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Histone, Variant
from tools.browse_service import *

import json
import os

import logging
from datetime import date, datetime
from cProfile import Profile

class Command(BaseCommand):
    help = 'Creating histone types and variants for the seed sequences'

    # Logging info
    logging.basicConfig(filename=os.path.join(LOG_DIRECTORY, "buildvarianttypes.log"),
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
        self.log.info('===             buildvarianttypes START             ===')
        self.log.info('=======================================================')
        self.start_time=datetime.now()

        #Clean the DB, removing all types and variants
        Variant.objects.all().delete()
        Histone.objects.all().delete()

        # Populate our Histone types table add descriptions
        self.create_histone_types()

        # Populate our Histone variants table add descriptions
        self.create_histone_variants()

        self.log.info('=======================================================')
        self.log.info('===     buildvarianttypes SUCCESSFULLY finished     ===')
        self.log.info('=======================================================')

    def handle(self, *args, **options):
        if options['profile']:
            profiler = Profile()
            profiler.runcall(self._handle, *args, **options)
            profiler.print_stats()
        else:
            self._handle(*args, **options)

    def create_histone_types(self):
        """Create basic histone types"""
        with open(VARIANTS_JSON) as f:
            variants_list = json.loads(f.read())
        for i in variants_list['tree'].keys():
            if i=='H1': continue
            obj,created = Histone.objects.get_or_create(id=i,taxonomic_span="Eukaryotes",\
                      description="Core histone")
            if created:
                self.log.info("Histone {} type was created in database.".format(obj.id))

        obj,created = Histone.objects.get_or_create(id="H1",taxonomic_span="Eukaryotes",\
                      description="Linker histone")
        if created:
            self.log.info("Histone {} type was created in database.".format(obj.id))

    # def get_variants(self):
    #     """Get iterator for variants with its histone type"""
    #     with open(VARIANTS_JSON) as f:
    #         variants_list = json.loads(f.read())
    #     for histtype in variants_list['tree'].keys():
    #         for var in variants_list['tree'][histtype].keys():
    #             yield histtype, var

    def create_histone_variants(self):
        """Create variants (including generics for each histone type) listed in variants_list.json"""
        with open(VARIANTS_JSON) as f:
            variants_classification = json.loads(f.read())
        variants_list = variants_classification['tree']
        for hist_type_pos in variants_list.keys():
            for variant_name in variants_list[hist_type_pos].keys():
                variant_model, create = Variant.objects.get_or_create(id=variant_name, hist_type_id=hist_type_pos)
                if create:
                    self.log.info("Created {} variant model in database".format(variant_model.id))
                for subvar_name in variants_list[hist_type_pos][variant_name].keys():
                    variant_model, create = Variant.objects.get_or_create(id=subvar_name, hist_type_id=hist_type_pos, parent_id=variant_name)
                    if create:
                        self.log.info("Created {} subvariant model in database".format(variant_model.id))
