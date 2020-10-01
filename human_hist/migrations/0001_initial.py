# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('browse', '__first__'),
    ]

    operations = [
        migrations.CreateModel(
            name='Histone_Human_cancer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('case_id', models.CharField(max_length=100)),
                ('sequencing_center', models.CharField(max_length=100)),
                ('cancer', models.CharField(max_length=100)),
            ],
        ),
        migrations.CreateModel(
            name='Histone_Human_genes',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('hgnc_symbol', models.CharField(max_length=25)),
                ('prev_hgnc_symb', models.CharField(max_length=25)),
                ('ncbi_gene_id', models.IntegerField()),
                ('ensg', models.CharField(max_length=100)),
                ('expr_timing', models.CharField(max_length=25)),
                ('expr_pattern', models.CharField(max_length=25)),
                ('biotype', models.CharField(max_length=25)),
                ('bona_fidecanonical', models.CharField(max_length=25)),
                ('pmids', models.CommaSeparatedIntegerField(max_length=250)),
                ('variant', models.ForeignKey(related_name='human_genes', to='browse.Variant')),
            ],
        ),
        migrations.CreateModel(
            name='Histone_Human_mutations',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('aa_change', models.CharField(max_length=100)),
                ('mutation_type', models.CharField(max_length=100)),
                ('ref_allele', models.CharField(max_length=100)),
                ('var_allele', models.CharField(max_length=100)),
                ('var_allele_freq', models.FloatField()),
                ('case', models.ManyToManyField(related_name='human_mutations', to='human_hist.Histone_Human_cancer')),
                ('gene', models.ManyToManyField(related_name='human_mutations', to='human_hist.Histone_Human_genes')),
            ],
        ),
        migrations.CreateModel(
            name='Histone_Human_proteins',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('enst', models.CharField(max_length=100)),
                ('refseq_transcript_id', models.CharField(max_length=100)),
                ('refseq_protein_id', models.CharField(max_length=100)),
                ('prot_lenght', models.IntegerField()),
                ('isoform', models.CharField(max_length=40)),
                ('gene', models.ManyToManyField(related_name='human_proteins', to='human_hist.Histone_Human_genes')),
            ],
        ),
    ]
