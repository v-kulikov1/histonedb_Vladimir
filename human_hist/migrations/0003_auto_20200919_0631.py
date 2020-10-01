# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('human_hist', '0002_auto_20200918_1216'),
    ]

    operations = [
        migrations.AlterField(
            model_name='histone_human_mutations',
            name='ref_allele',
            field=models.CharField(max_length=250),
        ),
        migrations.AlterField(
            model_name='histone_human_mutations',
            name='var_allele',
            field=models.CharField(max_length=150),
        ),
    ]
