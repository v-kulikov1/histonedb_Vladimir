# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('human_hist', '0001_initial'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='Histone_Human_cancer',
            new_name='Histone_Human_cancers',
        ),
    ]
