import math

from skmultiflow.drift_detection.base_drift_detector import BaseDriftDetector

class ANGLE(BaseDriftDetector):
    
	def __init__(self, delta = 0.05, block_size = 32, epsilon_prime = 0.01, alpha = 0.8, term = 75):
		self._delta = delta * 1.0
		self._block_size = block_size
		self.window = AngleWindow(block_size, 1, 1, epsilon_prime, alpha, term)
		self._element_count = 0
		self._drift_location = 0
		self._drift_angle = 0.0

	def detected_change(self):
		if self._element_count % self._block_size == 0 and self.window.get_block_count() >= 2:
			bln_reduce_width = True

			while bln_reduce_width:
				bln_reduce_width = False

				n1 = 0
				n0 = self.window.get_width()
				u1 = 0.0
				u0 = self.window.get_total()

				# new add
				block_means = [0.0] * self.window.get_block_count()
				block_sizes = [0] * self.window.get_block_count()
				cursor_index = self.window.get_block_count() - 1
				drift_location_index = cursor_index

				cursor = self.window.get_tail()

				while cursor.get_previous():
					#new
					block_means[cursor_index] = cursor.get_mean()
					block_sizes[cursor_index] = cursor.get_item_count()
					cursor_index -= 1

					n0 -= cursor.get_item_count()
					n1 += cursor.get_item_count()
					u0 -= cursor.get_total()
					u1 += cursor.get_total()

					diff = abs(u1 / n1 - u0 / n0)

					if diff > self.__get_ADWIN_bound(n0, n1):
						#new
						drift_location_index = cursor_index

						bln_reduce_width = True
						self.window.set_head(cursor)
						while cursor.get_previous():
							cursor = cursor.get_previous()
							#new
							block_means[cursor_index] = cursor.get_mean()
							block_sizes[cursor_index] = cursor.get_item_count()
							n1 += cursor.get_item_count()
							cursor_index -= 1

							self.window.set_width(self.window.get_width() - cursor.get_item_count())
							self.window.set_total(self.window.get_total() - cursor.get_total())
							self.window.set_variance(self.window.get_variance() - cursor.get_variance())
							self.window.set_block_count(self.window.get_block_count() - 1)

						self.window.get_head().set_previous(None)

						#new
						i_x = 0.0
						i_y = block_means[0]
						k_x = n1
						k_y = block_means[drift_location_index]
						j_x = 0.0
						j_y = i_y
						

						min_angle = 360.0
						angles = Buffer(2)
						angles_hist = [0.0] * (len(block_means)-1)

						for j in range(1, drift_location_index):
							j_x = j_x + block_sizes[j-1]
							j_y = block_means[j]

							angle = self.__get_angle(i_x, i_y, j_x, j_y, k_x, k_y)
							angles.add(angle)
							angles_hist[j] = angle

							if angles.get_mean() <= min_angle:
								min_angle = angle
								self._drift_location = int(k_x - j_x)
								self._drift_angle = min_angle

						return True

					cursor = cursor.get_previous()
		return False
	
	#new
	def __get_angle(self, i_x, i_y, j_x, j_y, k_x, k_y):

		u_x = j_x - i_x
		u_y = j_y - i_y
		v_x = j_x - k_x
		v_y = j_y - k_y

		u_mag = math.sqrt(u_x * u_x + u_y * u_y)
		v_mag = math.sqrt(v_x * v_x + v_y * v_y)

		u_dot_v = u_x * v_x + u_y * v_y

		cos_angle = u_dot_v / (u_mag * v_mag)
		if cos_angle > 1:
			cos_angle = 1.0
		elif cos_angle < -1:
			cos_angle = -1.0
		angle = math.acos(cos_angle)

		return angle * (180 / math.pi)


	def __get_ADWIN_bound(self, n0, n1):

		n0 = n0 * 1.0
		n1 = n1 * 1.0
		n = n0 + n1
		dd = math.log(2 * math.log(n) / self._delta)
		v = self.window.get_variance() / self.window.get_width()
		m = 1 / n0 + 1 / n1
		epsilon = math.sqrt(2 * m * v * dd) + 2 / 3 * dd * m

		return epsilon

	def add_element(self, value):
		self.window.add_transaction(value)
		self._element_count += 1

	def get_drift_location(self):
		return self._drift_location

	def get_drift_angle(self):
		return self._drift_angle



class AngleWindow(object):

	""" A linked-list like SeedWindow object for SEED algorithm.
    
    Used for storing SEED's block list. Is composed of SeedBlock objects. 
    Acts as a linked list, where each element points to its predecessor 
    and successor.
    """

	def __init__(self, block_size = None, decay_mode = None, compression_mode = None, epsilon_prime = None, alpha = None, compression_term = None):
		self.clear()
		self._block_size = block_size
		self._decay_mode = decay_mode if decay_mode else 1
		self._compression_mode = compression_mode if compression_mode else 1
		self._epsilon_prime = epsilon_prime * 1.0 if epsilon_prime else 0.0
		self._alpha = alpha * 1.0 if alpha else 0.0
		self._decay_compression_count = 0
		self._linear_fixed_term_size = compression_term if compression_term else 50 
		self.add_block_to_head(AngleBlock(block_size))

	def clear(self):
		self._head = None
		self._tail = None
		self._width = 0
		self._block_count = 0
		self._total = 0.0
		self._variance = 0.0

	def add_transaction(self, value):

		value = value * 1.0

		if self._tail.is_full():
			if self._compression_mode == 1:
				if self._tail.get_previous() and self._decay_compression_count > self._linear_fixed_term_size:
					self._decay_compression_count = 0
					cursor = self._tail
					epsilon = 0.0
					i = 0

					while cursor and cursor.get_previous():
						n0 = cursor.get_item_count() * 1.0
						n1 = cursor.get_previous().get_item_count() * 1.0
						u0 = cursor.get_total() * 1.0
						u1 = cursor.get_previous().get_total() * 1.0

						diff = abs(u1 / n1 - u0 / n0)

						if self._decay_mode == 1:
							epsilon += self._epsilon_prime * self._alpha
						elif self._decay_mode == 2:
							epsilon = self._epsilon_prime * math.pow(1 + self._alpha, i)

						if diff < epsilon:
							self.compress_block(cursor)

						cursor = cursor.get_previous()

						i += 1

			self.add_block_to_tail(AngleBlock(self._block_size))
			self._decay_compression_count += 1

		self._tail.add(value)
		self._total += value
		self._width += 1

		if self._width > 2:
			inc_variance = (self._width - 1) * (value - self._total / (self._width - 1)) * (value - self._total / (self._width - 1)) / self._width
			self._variance += inc_variance
			self._tail.set_variance(self._tail.get_variance() + inc_variance)

	def compress_block(self, block):
		block.get_previous().set_total(block.get_total() + block.get_previous().get_total())
		block.get_previous().set_item_count(block.get_item_count() + block.get_previous().get_item_count())
		block.get_previous().set_variance(block.get_variance() + block.get_previous().get_variance())
		block.get_previous().set_block_size(block.get_block_size() + block.get_previous().get_block_size())

		if block.get_next():
			block.get_previous().set_next(block.get_next())
			block.get_next().set_previous(block.get_previous())
		else:
			block.get_previous().set_next(None)
			self._tail = block.get_previous()

		self._block_count -= 1

	#def check_homogeneity(block):

	#	diff = abs(block.get_mean() - block.get_previous.get_mean())
	#	epsilon_prime = self.__get_ADWIN_bound(block.get_item_count(), block.get_previous().get_item_count())
	#	return diff < epsilon_prime

	#def __get_ADWIN_bound(self, n0, n1):

	#	n0 = n0 * 1.0
	#	n1 = n1 * 1.0
	#	n = n0 + n1
	#	dd = math.log(2 * math.log(n) / 0.99)
	#	v = self._variance / self._width
	#	m = 1 / n0 + 1 / n1
	#	epsilon = math.sqrt(2 * m * v * dd) + 2 / 3 * dd * m

	#	return epsilon


	def add_block_to_head(self, block):

		if self._head:
			block.set_next(self._head)
			self._head.set_previous(block)
			self._head = block
		else:
			self._head = block
			self._tail = block

		self._block_count += 1

	def remove_block(self, block):
		self._width -= block.get_item_count()
		self._total -= block.get_total()
		self._variance -= block.get_variance()
		self._block_count -= 1

		if block.get_previous() and block.get_next():
			block.get_previous().set_next(block.get_next())
			block.get_next().set_previous(block.get_previous())
			block.set_next(None)
			block.set_previous(None)
		elif not block.get_previous() and block.get_next():
			block.get_next().set_previous(None)
			self._head = block.get_next()
			block.set_next(None)
		elif block.get_previous() and not block.get_next():
			block.get_previous().set_next(None)
			self._tail = block.get_previous()
			block.set_previous(None)
		else:
			self._head = None
			self._tail = None

	def add_block_to_tail(self, block):

		if self._tail:
			block.set_previous(self._tail)
			self._tail.set_next(block)
			self._tail = block
		else:
			self._head = block
			self._tail = block

		self._block_count += 1


	def set_block_count(self, value):
		self._block_count = value

	def get_block_count(self):
		return self._block_count

	def set_width(self, value):
		self._width = value

	def get_width(self):
		return self._width

	def set_head(self, head):
		self._head = head

	def get_head(self):
		return self._head

	def set_tail(self, tail):
		self._tail = tail

	def get_tail(self):
		return self._tail

	def set_total(self, value):
		self._total = value

	def get_total(self):
		return self._total

	def set_variance(self, value):
		self._variance = value

	def get_variance(self):
		return self._variance

	def set_block_size(self, value):
		self._block_size = value if value > 32 else 32

	def get_block_size(self):
		return self._block_size

	def set_epsilon_prime(self, value):
		self._epsilon_prime = value

	def get_epsilon_prime(self):
		return self._epsilon_prime

	def set_alpha(self, value):
		self._alpha = value

	# new add
	def get_alpha(self):
		return self._alpha

	def set_compression_term(self, value):
		self._linear_fixed_term_size = value

	# new add
	def get_compression_term(self):
		return self._linear_fixed_term_size

class AngleBlock(object):
	""" SeedBlock to be used by the SeedWindow object.
    
    The SeedBlock object, alongside the SeedWindow object, are the two main data 
    structures used for storing the relevant statistics for the SEED
    algorithm for change detection.
    
    Parameters
    ----------
    block_size: integer number
        Block size for the new created block
    block: SeedBlock object
        Reference to the current SeedBlock in the BlockWindow
    
    """

	def __init__(self, block_size = None, block = None):
		self._next = None if block_size else block.get_next()
		self._previous = None if block_size else block.get_previous()
		self._block_size = block_size if block_size else block._block_size
		self._total = 0 if block_size else block._total
		self._variance = 0 if block_size else block._variance
		self._item_count = 0 if block_size else block._item_count

	def set_next(self, block):
		self._next = block

	def get_next(self):
		return self._next

	def set_previous(self, block):
		self._previous = block

	def get_previous(self):
		return self._previous

	def set_block_size(self, block_size):
		self._block_size = block_size

	def get_block_size(self):
		return self._block_size

	def set_total(self, value):
		self._total = value

	def get_total(self):
		return self._total

	def set_item_count(self, value):
		self._item_count = value

	def get_item_count(self):
		return self._item_count

	def set_variance(self, value):
		self._variance = value

	def get_variance(self):
		return self._variance

	def add(self, value):
		self._item_count += 1
		self._total += value

	def is_full(self):
		return self._item_count == self._block_size

	def get_mean(self):
		return 1.0 * self._total / self._item_count  # exception

class Buffer(object):

	def __init__(self, size):

		self._buffer = [0.0] * size
		self._size = size
		self._sliding_index = 0
		self._is_full = False
		self._total = 0.0

	def add(self, value):

		if self._sliding_index == self._size:
			self._is_full = True
			self._sliding_index = 0

		removed = self._buffer[self._sliding_index]
		self._total -= removed

		self._buffer[self._sliding_index] = value
		self._sliding_index += 1
		self._total += value

		return removed if self._is_full else -1

	def get_buffer(self):

		return self._buffer

	def get_mean(self):

		return self._total / self._size if self._is_full else self._total / self._sliding_index

	def get_std_dev(self):

		return calculateStdev(self._buffer, self.get_mean())

	def __calculateStdev(self, times, mean):

		sum_all = 0.0
		count = 0

		for time in times:
			if time > 0.0:
				count += 1
				sum_all += math.pow(time-mean, 2)

		return math.sqrt(sum_all / count) 


	def is_full(self):

		return self._is_full

	def clear(self):
		self._buffer = [0.0] * size
		self._sliding_index = 0
		self._is_full = False
		self._total = 0.0















