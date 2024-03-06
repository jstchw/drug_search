import { motion, Variants } from 'framer-motion';
import { ControlledMenu } from '@szhsin/react-menu';
import React, { ReactNode } from 'react';

const OFFSET_X = 50;
const OFFSET_Y = 50;

const menuVariants: Variants = {
  open: {
    opacity: 1,
    scale: 1,
    transition: { type: 'spring', stiffness: 260, damping: 20, delay: 0.1 },
    x: 0,
    y: 0,
  },
  closed: { opacity: 0, scale: 0.9, x: OFFSET_X, y: OFFSET_Y },
};

interface AnimatedMenuProps {
  isOpen: boolean;
  anchorPoint: { x: number; y: number };
  onClose: () => void;
  children: ReactNode;
}

const AnimatedMenu: React.FC<AnimatedMenuProps> = ({ isOpen, anchorPoint, onClose, children }) => {
  const adjustedAnchorPoint = { x: anchorPoint.x + OFFSET_X, y: anchorPoint.y + OFFSET_Y };

  return (
    <motion.div
      variants={menuVariants}
      initial="closed"
      animate={isOpen ? 'open' : 'closed'}
      style={{ position: 'absolute', zIndex: 1000 }}
    >
      <ControlledMenu anchorPoint={adjustedAnchorPoint} state={isOpen ? 'open' : 'closed'} onClose={onClose}>
        {children}
      </ControlledMenu>
    </motion.div>
  );
};

export default AnimatedMenu;
