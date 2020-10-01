from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^human2/$', views.human_histones),
]