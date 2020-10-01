# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('human_hist', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='histone_human_mutations',
            name='case',
            field=models.ManyToManyField(related_name='human_mutations_2', to='human_hist.Histone_Human_cancers'),
        ),
    ]
