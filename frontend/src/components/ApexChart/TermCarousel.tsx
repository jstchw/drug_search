import { ApexOptions } from 'apexcharts';
import { URLParams } from '../../types';
import ReactApexChart from 'react-apexcharts';
import useGeneralOptionsStore from '../../stores/generalOptionsStore';
import React from 'react';
import ReactWordcloud, { OptionsProp, CallbacksProp } from 'react-wordcloud';
import { Carousel } from 'react-responsive-carousel';
import { Nav, OverlayTrigger, Popover, Pagination, Row, Col } from 'react-bootstrap';
import { Cloud, List, Pill } from '@phosphor-icons/react';
import { SealWarning, ChartLine, SmileyNervous } from '@phosphor-icons/react';
import { AnimatePresence, motion } from 'framer-motion';
import useDemographicStore from '../../stores/demographicStore';
import { useUrlParams } from '../../hooks/useUrlParams';
import { useTermData } from '../../hooks/useTermData';
import { searchTypes } from '../../constants';
import { capitalizeFirstLetter, getColorFromPercentage, valueToPercentage } from '../../utils/utils';
import _ from 'lodash';
import 'tippy.js/dist/tippy.css';
import 'tippy.js/animations/scale.css';
import { MenuItem } from '@szhsin/react-menu';
import ContextMenu from '../ContextMenu/ContextMenu';
import '@szhsin/react-menu/dist/index.css';
import useArticleStore from '../../stores/articleStore';

const cloudOptions: OptionsProp = {
  enableTooltip: true,
  enableOptimizations: true,
  deterministic: true,
  fontFamily: 'trebuchet ms',
  fontSizes: [10, 60],
  fontStyle: 'normal',
  fontWeight: 'normal',
  padding: 1,
  rotations: 0,
  scale: 'log',
  spiral: 'archimedean',
  transitionDuration: 1000,
};

export const getChartWarning = (params: URLParams) => {
  const popover = (
    <Popover>
      <Popover.Body>
        <div className={'d-flex align-items-center justify-content-center'}>
          <ChartLine weight={'light'} />
          <div className={'vr mx-2'} />
          Correlation does not imply causation.
        </div>
      </Popover.Body>
    </Popover>
  );

  return (
    <OverlayTrigger trigger={['hover', 'focus']} placement={'bottom'} overlay={popover}>
      <div style={{ cursor: 'default' }}>
        {params.searchBy !== 'side_effect' ? (
          <span className={'d-inline-flex align-items-center'}>
            <SmileyNervous className={'text-secondary'} weight={'light'} />
            <SealWarning className={'text-secondary'} weight={'light'} />
            <div className={'vr mx-2'} />
            <span>Side Effects</span>
          </span>
        ) : (
          <span className={'d-inline-flex align-items-center'}>
            <Pill className={'text-secondary'} weight={'light'} />
            <SealWarning className={'text-secondary'} weight={'light'} />
            <div className={'vr mx-2'} />
            <span>Substances</span>
          </span>
        )}
      </div>
    </OverlayTrigger>
  );
};

// Determine the display type based on the searchBy parameter
const determineDisplayType = (searchBy: string) => {
  if (searchTypes[0] !== undefined && searchTypes[2] !== undefined) {
    // If the search parameter is not a side effect (generic or brand name)
    // Return the opposite of the search parameter (to query API for side effects)
    // Otherwise, return the generic name (to query API for generic names)
    if (searchBy !== searchTypes[2].param) {
      return searchTypes[2].param;
    } else {
      return searchTypes[0].param;
    }
  } else {
    return '';
  }
};

interface TermCarouselProps {
  noFilterRequest?: boolean;
  onRender: () => void;
  source: string;
}

const TermCarousel: React.FC<TermCarouselProps> = ({ noFilterRequest = false, onRender, source }) => {
  const {
    params: { searchBy },
  } = useUrlParams();

  const theme = useGeneralOptionsStore((state) => state.theme);
  const [carouselIndex, setCarouselIndex] = React.useState<number>(0);

  const { data, isError } = useTermData(source, noFilterRequest);

  const setShowDemographic = useDemographicStore((state) => state.setShowDemographic);
  const setDemographicTerm = useDemographicStore((state) => state.setDemographicTerm);
  const setDemographicType = useDemographicStore((state) => state.setDemographicType);

  /* Context menu declarations */
  const [isContextMenuOpen, setIsContextMenuOpen] = React.useState<boolean>(false);
  const [contextMenuAnchor, setContextMenuAnchor] = React.useState({ x: 0, y: 0 });
  const [contextMenuSelectedTerm, setContextMenuSelectedTerm] = React.useState<string>('');

  const openDemographicModal = () => {
    setDemographicTerm(contextMenuSelectedTerm);
    setDemographicType(determineDisplayType(searchBy));
    setShowDemographic(true);
  };

  const addArticleTerm = useArticleStore((state) => state.addArticleTerm);

  // Callback to parent component to indicate that the component has rendered
  React.useEffect(() => {
    if (data && !isError) {
      onRender();
    }
  }, [data, isError]);

  const [currentChartPage, setCurrentChartPage] = React.useState<number>(0);
  const itemsPerPage = 10;

  const handleChartPageChange = (pageNumber: number) => {
    setCurrentChartPage(pageNumber);
  };

  if (!data || isError) {
    return null;
  }

  // Pagination items
  const seriesForCurrentPage = data?.series.map((series) => ({
    ...series,
    data: series.data.slice(currentChartPage * itemsPerPage, (currentChartPage + 1) * itemsPerPage),
  }));
  const categoriesForCurrentPage = data?.categories.slice(
    currentChartPage * itemsPerPage,
    (currentChartPage + 1) * itemsPerPage
  );

  const totalPageCount = Math.ceil(data.categories.length / itemsPerPage);
  const paginationItems = [];
  for (let number = 0; number < totalPageCount; number++) {
    paginationItems.push(
      <Pagination.Item key={number} active={number === currentChartPage} onClick={() => handleChartPageChange(number)}>
        {number + 1}
      </Pagination.Item>
    );
  }

  /* Apex Chart declarations */
  const chartOptions: ApexOptions = {
    colors: [
      function ({ value }: { value: number }) {
        const percentage = valueToPercentage(value, data.total_count);
        return getColorFromPercentage(percentage);
      },
    ],
    theme: {
      mode: theme,
    },
    legend: {
      show: false,
    },
    chart: {
      width: '100%',
      toolbar: {
        show: false,
        tools: {
          zoom: false,
          zoomin: false,
          zoomout: false,
        },
      },
      background: theme === 'dark' ? '#212529' : '',
      dropShadow: {
        enabled: true,
        top: 1,
        left: 1,
        blur: 2,
        color: theme === 'dark' ? '#000' : '#000',
        opacity: 0.2,
      },
      events: {
        dataPointSelection: (event, _, config) => {
          setContextMenuSelectedTerm(capitalizeFirstLetter(config.w.config.labels[config.dataPointIndex]));
          setContextMenuAnchor({ x: event.clientX, y: event.clientY });
          setIsContextMenuOpen(true);
        },
      },
    },
    states: {
      active: {
        filter: {
          type: 'none',
        },
      },
    },
    plotOptions: {
      bar: {
        horizontal: true,
        distributed: true,
        barHeight: '80%',
        borderRadius: 4,
        borderRadiusApplication: 'end',
      },
    },
    grid: {
      show: false,
    },
    labels: categoriesForCurrentPage,
    yaxis: {
      labels: {
        show: true,
        style: {
          fontSize: '14px',
          colors: theme === 'dark' ? '#ACB5BD' : '',
        },
        formatter: (value: number) => {
          const stringValue = String(value);
          return stringValue.charAt(0).toUpperCase() + stringValue.slice(1).toLowerCase();
        },
      },
    },
    xaxis: {
      labels: {
        show: false,
      },
    },
    tooltip: {
      enabled: isContextMenuOpen ? false : true,
      y: {
        formatter: (value: number) => {
          return `${valueToPercentage(value, data.total_count).toPrecision(3)}% from ${data.total_count.toLocaleString()}`;
        },
      },
    },
    dataLabels: {
      enabled: true,
      // Converts the value to a percentage
      // formatter: (value: number) => {
      //   return `${((value / data.total_count) * 100).toPrecision(3)}%`;
      // },
      formatter: (value: number) => {
        return `${value.toLocaleString()}`;
      },
    },
  };

  /* Word Cloud declarations*/
  const cloudData = data.categories.map((category, index) => {
    const value = data.series[0]?.data[index] as number;
    return {
      text: category.toUpperCase(),
      value: value,
    };
  });

  const cloudCallbacks: CallbacksProp = {
    getWordColor: (word: { text: string; value: number }) => {
      const maxFrequency = Math.max(...(cloudData.map((w) => w.value) as number[]));
      const frequencyRatio = word.value / maxFrequency;

      if (theme === 'dark') {
        const darkGray = [172, 181, 189];
        const brightRed = [255, 0, 0];

        const interpolatedColor = darkGray.map((start, i) => {
          const end = brightRed[i];
          if (end !== undefined) {
            return Math.floor(start + frequencyRatio * (end - start));
          } else {
            return start;
          }
        });

        return `rgb(${interpolatedColor.join(', ')})`;
      } else {
        const redComponent = Math.floor(frequencyRatio * 255);
        return `rgb(${redComponent}, 0, 0)`;
      }
    },
    onWordClick: (word: { text: string }, event) => {
      setContextMenuSelectedTerm(capitalizeFirstLetter(word.text));
      setContextMenuAnchor({ x: event!.clientX, y: event!.clientY });
      setIsContextMenuOpen(true);
    },
  };

  return (
    <>
      <AnimatePresence>
        <ContextMenu
          anchorPoint={contextMenuAnchor}
          isOpen={isContextMenuOpen}
          onClose={() => setIsContextMenuOpen(false)}
        >
          <MenuItem onClick={openDemographicModal} className={(state) => (state.hover ? 'text-black' : '')}>
            Demographic breakdown
          </MenuItem>
          <MenuItem
            onClick={() => addArticleTerm(contextMenuSelectedTerm, determineDisplayType(searchBy), true)}
            className={(state) => (state.hover ? 'text-black' : '')}
          >
            Add to article filters
          </MenuItem>
        </ContextMenu>
      </AnimatePresence>
      <motion.div
        layout
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        transition={{
          type: 'spring',
          stiffness: 260,
          damping: 20,
        }}
        exit={{ opacity: 0 }}
      >
        <Row className={'mt-1'}>
          <Col className={'mb-3 text-center'}>
            <span className={'text-secondary'}>
              {data.total_count.toLocaleString()} events collected for {data.categories.length} terms
            </span>
          </Col>
        </Row>
        <Nav variant="tabs" defaultActiveKey={carouselIndex} className={'mt-3 z-index-n-1'}>
          <Nav.Item>
            <Nav.Link className={'d-flex align-items-center'} eventKey="0" onClick={() => setCarouselIndex(0)}>
              <List weight={'light'} />
              <div className={'vr mx-2'} />
              Term Chart
            </Nav.Link>
          </Nav.Item>
          <Nav.Item>
            <Nav.Link className={'d-flex align-items-center'} eventKey="1" onClick={() => setCarouselIndex(1)}>
              <Cloud weight={'light'} />
              <div className={'vr mx-2'} />
              Word Cloud
            </Nav.Link>
          </Nav.Item>
        </Nav>
        <Carousel
          showThumbs={false}
          showIndicators={false}
          showArrows={false}
          showStatus={false}
          selectedItem={carouselIndex}
          swipeable={false}
        >
          <div>
            <ReactApexChart options={chartOptions} series={seriesForCurrentPage} type={'bar'} />
            <Pagination size={'lg'} className={'d-flex justify-content-center'}>
              {paginationItems}
            </Pagination>
          </div>
          <ReactWordcloud words={cloudData} options={cloudOptions} callbacks={cloudCallbacks} />
        </Carousel>
      </motion.div>
    </>
  );
};

export default TermCarousel;
