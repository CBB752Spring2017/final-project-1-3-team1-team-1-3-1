from django.conf.urls import url, include
from django.views.generic import ListView, DetailView
from network_analysis.models import Node

urlpatterns = [ url(r'^$', ListView.as_view(queryset=Node.objects.all().order_by("-degree_centrality"), template_name="network_analysis/network_analysis.html")),
               url(r'^(?P<pk>\d+$)',DetailView.as_view(model= Node, template_name = "network_analysis/detail.html"))
               ]
