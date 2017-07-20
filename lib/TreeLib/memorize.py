from collections import Hashable as hashable
from functools import partial

class memorize(object):
	"""Cache function output when it's called and return it
	later when the same function is called with the same input,
	in this case, memorize use a hash to determine value to reevalute
	"""
	def __init__(self, function):
		self.function=function
		self.cache={}

	def __call__(self, hashrep, *args, **kwargs):
		"""Call to memorize, (as decorator)"""

		if hashrep in self.cache:
			print "FUCK YEAHHHH"
			return self.cache[hashrep]

		elif not isinstance(hashrep, hashable) or hashrep is None:
			#hashrep is None or uncachable
			return self.function(*args, **kwargs)

		else:
			output = self.function(*args, **kwargs)
			self.cache[hashrep]=output
			return output

	def __repr__(self):
		"""Return cached data"""
		for hashrep, data in self.cache:
			print hashrep, '=============>\n', data

	def __get__(self, obj, objtype):
		"""Instance methods support"""
		return partial(self.__call__,obj )
