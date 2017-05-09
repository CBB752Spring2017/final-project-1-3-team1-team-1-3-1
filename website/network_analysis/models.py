from django.db import models

# Create your models here.

class Node(models.Model):
    name = models.CharField(max_length = 140)
    uniprot = models.CharField(max_length = 140)
    degree_centrality = models.DecimalField(max_digits = 50, decimal_places=10)
    betweenness_centrality = models.DecimalField(max_digits = 50, decimal_places=10)

    def __str__(self):
        return self.name
