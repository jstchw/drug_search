import { Nav, Col } from 'react-bootstrap';
import SearchBox from '../SearchBox/SearchBox';
import Header from '../Header/Header';
import { motion, useInView, AnimatePresence } from 'framer-motion';
import { FC, useState, useEffect, RefObject } from 'react';
import { isMobile } from 'react-device-detect';

interface ScrollbarProps {
  searchBoxSectionRef: RefObject<HTMLDivElement>;
}

const Scrollbar: FC<ScrollbarProps> = ({ searchBoxSectionRef }) => {
  const [showScrollbar, setShowScrollbar] = useState<boolean>(false);
  const isSearchBoxInView = useInView(searchBoxSectionRef);

  useEffect(() => {
    setShowScrollbar(!!!isSearchBoxInView);
  }, [isSearchBoxInView]);

  return (
    <>
    <AnimatePresence>
      {showScrollbar ? (
          <motion.div
            initial={{ y: -200 }}
            animate={{ y: 0 }}
            exit={{ y: -200 }}
            transition={{
              type: 'spring',
              stiffness: 300,
              damping: 24,
              duration: 0.8,
              delay: 0.1,
            }}
            className="fixed-top py-4 rounded shadow-lg bg-body"
          >
            <Nav>
              <Col 
                className={
                  !isMobile ? 'ms-4' : 'justify-content-center align-items-center d-flex mb-3'
                }
              >
                <Header />
              </Col>

              <Col className="d-flex justify-content-center">
                <SearchBox />
              </Col>

              <Col xs={4}></Col>
            </Nav>
          </motion.div>
      ) : undefined}
    </AnimatePresence>
    </>
  );
};

export default Scrollbar;