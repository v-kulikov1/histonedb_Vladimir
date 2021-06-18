from django.conf.urls import url
from django.templatetags.static import static
from . import views

app_name = 'human_hist'
urlpatterns = [
    url(r'^$', views.human_histones),
    url(r'hist_mutations', views.hist_mutations),
]